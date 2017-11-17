#!/usr/bin/env perl
# purpose: annotate HIV consensus sequence with gene/CDS regions
# report information about mutations with the reference and protein productivity

use Bio::Perl;
use Bio::DB::GenBank;
use Bio::SeqFeature::Generic;
use Bio::Location::Split;
use Bio::Tools::Run::Alignment::Clustalw;
use POSIX;
use Text::Table;
use Data::Printer;
use feature qw(say);

my @gene_names = ('gag','pol','env','envelope','vif','vpr','tat','rev','vpu','nef');

my $write_out, $outfilename;
if ( $ARGV[2] eq 'consensus' ) {
    $write_out = 0;
    $outfilename = ( split /.fasta$/, $ARGV[0])[0];
} else {
    $write_out = 1;
    $outfilename = ( split /_out_\S*.unpadded.fasta/, $ARGV[0])[0];
}

# PARAMETERS TO CHANGE:
my $min_length = 0.6; # minimum fraction of full protein length to be considered productive
my $max_gap = 10; # maximum gap size to be considered deletion rather than missing sequence (bp)

# read in consensus sequence
die "ERROR: Input consensus fasta ('$ARGV[0]') is not a valid file\n$!" if ! -e -s $ARGV[0];
my $sequence_id = ( split /\//, $outfilename )[-1];
my $consensus_file = Bio::SeqIO->new( -file => $ARGV[0],
                                      -format => "fasta" );
my $consensus = $consensus_file->next_seq();
#my $consensus_id = substr $consensus->id, 0, 30; # truncate to 30 characters to match clustalw
my $consensus_id = $consensus->id;
( my $clean_seq = uc($consensus->seq) ); #=~ tr/[NX]//d;  # remove N's (x's) from consensus
$consensus->seq($clean_seq);

# get genbank entry for best reference
say "consensus sequence id: " . $consensus->id;
my $ref = (split /_/, $consensus->id)[-2];
say "Genbank reference: " . $ref;
#my $ref = 'U51190';
my $io = Bio::SeqIO->new( -file => $ARGV[1] . "/" . $ref . ".gb",
                          -format => "genbank" );
my $ref_obj = $io->next_seq();
my $ref_seq = Bio::Seq->new( -id => $ref,
                             -seq => $ref_obj->seq );


# create genbank for concatenated sequence
my $concat = Bio::Seq->new( -id => $consensus->id,
                            -seq=> "" );

say 'creating mutation annotation and GenBank file for ' . $sequence_id;
open( my $outfile, '>', $outfilename.'_annotation.txt') or die "Could not open file '$outfilename'._annotation.txt $!";
say $outfile  '====================================================================================';
say $outfile  "SEQ_ID:\t" . $consensus->id;
say $outfile  "REF_ID:\t" . $ref_seq->id;

open( my $outseqfile, '>', $outfilename.'_seq.txt') or die "Could not open file '$outfilename'._seq.txt $!";
say $outseqfile  ">" . $consensus->id;

sub align {
    # align consensus sequence to best reference using clustalw (suppressing output)
    my ($ref_seq_obj, $consensus_seq_obj, $write_out) = @_;
    @params = ('quiet' => 1);
    $factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);

    my @sequences = ($ref_seq_obj, $consensus_seq_obj);
    my $aln = $factory->align(\@sequences);

    # output alignment
    if ($write_out) {
        $out = Bio::AlignIO->new( -file => ">analysis/".$sequence_id."/full_alignment.aln",
                                  -format => 'clustalw');
        $out->write_aln($aln);
    }
    return $aln;
}


sub check_frame {
    my ($cds_seq) = @_;
    my (@inframe_idx, @outframe_idx, @missing_idx);

    while ($cds_seq =~ /-+/g) {
        if ( length($&) % 3 == 0 ) {
            push @inframe_idx, ($-[0]..$+[0]);
        } elsif ( length($&) % 3 == 1 ) {
            push @outframe_idx, ($-[0]);
            push @inframe_idx, ($-[0]+1..$+[0]);
        } elsif ( length($&) % 3 == 2 ) {
            push @outframe_idx, ($-[0], $-[0]+1);
            push @inframe_idx, ($-[0]+2..$+[0]);
        }

        # call gaps of more than max_gap bp missing sequence instead of deletions
        if ( length($&) > $max_gap ) {
            push @missing_idx, ($-[0]..$+[0]);
        }
    }
    return (\@inframe_idx, \@outframe_idx, \@missing_idx);
}


sub find_mutations {
    my ($cds_aln, $cds_start, $cds_end, $cons_protein, $ref_protein, $gene_name, $write_out) = @_;

    my %mutations = (
        mismatches => 0,
        insertions => 0,
        deletions => 0,
    );
    if ($write_out){
        $out = Bio::AlignIO->new( -file => ">analysis/".$sequence_id."/".$gene_name."_alignment.aln",
                                  -format => 'clustalw');
        $out->write_aln($cds_aln);
    }
    my $first_pos = $cds_aln->column_from_residue_number($consensus_id, $cds_start);
    my $last_pos = $cds_aln->column_from_residue_number($consensus_id, $cds_end);
    #say $outfile  join "\t", $first_pos . $last_pos . $cds_aln->length;

    # get reference and consensus sequences for CDS, setting gaps as '-'
    my $ref_seq = $cds_aln->get_seq_by_id($ref);
    (my $ref_string = $ref_seq->seq) =~ tr/\./-/;

    my $consensus_seq = $cds_aln->get_seq_by_id($consensus_id);
    (my $consensus_string = $consensus_seq->seq) =~ tr/\./-/;

    my $ref_start = $ref_seq->location_from_column(1)->start;
    my $cons_start = $consensus_seq->location_from_column($first_pos)->start;

    # characterize indels as inframe (multiple of 3) or outframe (causes frameshift)
    my ($in, $out, $miss) = check_frame($consensus_string);
    my @inframe_dels = @$in;
    my @outframe_dels = @$out;
    my @missing_idx = @$miss;

    my ($in, $out, $miss) = check_frame($ref_string);
    my @inframe_ins = @$in;
    my @outframe_ins = @$out;

    my @mutation_lines;         # hold mutation print lines for output later
    my $offset = 0;             # frame shift counter
    BASE: foreach my $pos ( $first_pos..$last_pos ) { # more precise indexing?
    #foreach my $pos ( 1..$cds_aln->length ) {
        # compare base by base
        my $ref_base = substr $ref_string, $pos-1, 1;
        my $cons_base = substr $consensus_string, $pos-1, 1;

        if ($ref_base ne $cons_base) {
            # get corresponding amino acid
            my $ref_res = $ref_seq->location_from_column($pos);
            my $cons_res = $consensus_seq->location_from_column($pos);
            my $ref_idx = $ref_res->start - $ref_start;
            my $cons_idx = $cons_res->start - $cons_start;
            #say $outfile  $ref_idx . "\t" . $cons_idx;

            my $ref_aa = ( $ref_res->location_type eq 'EXACT' ) ? substr $ref_protein, floor($ref_idx/3), 1 : '-';
            my $cons_aa = ( $cons_res->location_type eq 'EXACT' ) ? substr $cons_protein, floor($cons_idx/3), 1 : '-';

            # distinguish between deletion and missing sequence
            if ( grep { $_ eq $pos-1 } @missing_idx ){
                my @outline = ("\t\t", $pos, $ref_base, $cons_base, $ref_aa, "missing sequence");
                push @mutation_lines, \@outline;
                next BASE;
            }
            #if ($pos >= $first_pos && $pos <= $last_pos) { # distinguish between true gaps and missing sequence
            # determine type of mutation
            if ( $cons_base eq '-' ) { # deletion
                $mutations{deletions} += 1;
                if ( grep { $_ eq $pos-1 } @outframe_dels ) {
                    $offset -= 1;
                    $offset = 0 if $offset == -3;
                }
            } elsif ( $ref_base eq '-' ) { # insertion
                $mutations{insertions} += 1;
                if ( grep { $_ eq $pos-1 } @outframe_ins ) {
                    $offset += 1;
                    $offset = 0 if $offset == 3;
                }
            } else { # mismatch
                $mutations{mismatches} += 1;
            }

            # determine if nucleotide mutation is synonymous (* = synonymous)
            my $equiv = ( $cons_aa eq $ref_aa ) ? '*' : '';
            my $frameshift = ($offset != 0) ? $offset : '';
            my @outline = ("\t\t", $pos, $ref_base, $cons_base, $ref_aa, $cons_aa, $equiv, $frameshift);
            push @mutation_lines, \@outline;
        }
    }
    return (\%mutations, \@mutation_lines);
}


sub main {
    my $aln = align($ref_seq, $consensus, $write_out);
    my $aln = $aln->set_new_reference($ref);
    my %seqs;

    FEAT: for my $feat_object ($ref_obj->get_SeqFeatures) {
        my ($gene_name, $ref_protein, $cds_aln, $cds_string, $first_col, $last_col);
        my $splitlocation = Bio::Location::Split->new();
        my ($feat_start, $add_loc, $split_gene);

        # extract information about CDS regions
        if ($feat_object->primary_tag eq 'CDS') {
            #print "coords: ", $feat_object->start, "-", $feat_object->end, "\n";
            for my $tag ($feat_object->get_all_tags) {
                my @tag_vals = $feat_object->get_tag_values($tag);
                if ( ($tag eq 'note' || $tag eq 'product') && ! $gene_name ) {
                    $gn = lc ((split / /, $tag_vals[0])[0]);
                    if (grep { $_ eq $gn } @gene_names) {
                        $gene_name = $gn;
                    }
                } elsif ($tag eq 'gene') {
                    $gene_name = lc ($tag_vals[0]);
                } elsif ($tag eq 'translation') {
                    $ref_protein = $tag_vals[0]."*";
                }
            }
            # continue only if gene name and translation are found (feature is valid CDS)
            unless ( $ref_protein && $gene_name ) {
                next FEAT;
            }

            say $outfile '====================================================================================';
            say $outfile "Gene: ", $gene_name;
            #say  ' seq: ', $cds_protein;
            #say ' ref: ', $ref_protein;
            #say "Gene: ", $gene_name;

            # splice alignment to get current CDS region only
            if ($feat_object->location->isa('Bio::Location::SplitLocationI')) {     # spliced CDS
                say 'gene: ', $gene_name;
                SUBLOC: for my $location ($feat_object->location->sub_Location) {
                    $first_col = $aln->column_from_residue_number($ref, $location->start);
                    $last_col = $aln->column_from_residue_number($ref, $location->end);

                    $add_loc = 1 if $split_gene;
                    $split_gene = $gene_name;
                    $cds_aln = $aln->slice($first_col, $last_col, 0);
                    if ($cds_aln->num_sequences != 2) {
                        say $outfile "no sequence";
                        next FEAT;
                    }

                    # extract consensus CDS and reference CDS sequences
                    my $cds_seq = $cds_aln->get_seq_by_id($consensus_id);
                    ( my $tmp_seq = $cds_seq->seq ) =~ tr/\.//d;    # remove gaps from consensus
                    $cds_string = $cds_string . $tmp_seq;    # remove gaps from consensus

                    #say $gene_name . "\t\t\t" . "coords: ", $feat_object->start, "-", $feat_object->end;
                    my $orig_start = index($clean_seq, $tmp_seq) + 1;
                    my $orig_end = $orig_start + length($tmp_seq) - 1;

                    $splitlocation->add_sub_Location(Bio::Location::Simple->new( -start => $orig_start,
                                                                                 -end => $orig_end ) );

                    if ($add_loc == 1) {    # seen all 'exons' of CDS; re-align and create CDS feature for genbank file
                        my $cds_protein = translate_as_string($cds_string);
                        say 'cds nucleotides: ', $cds_string;
                        say 'ref nucleotides: ', $feat_object->spliced_seq->seq;
                        #say 'cds protein: ', $cds_protein;
                        #say 'ref protein: ', $ref_protein;

                        # report mutations
                        my $cds_ref_seq = Bio::Seq->new( -id => $ref,
                                                         -seq => $feat_object->spliced_seq->seq );
                        my $cds_cons_seq = Bio::Seq->new( -id => $consensus_id,
                                                         -seq => $cds_string );
                        #my $realn = align($cds_ref_seq, $cds_cons_seq, $gene_name);
                        my $realn = align($cds_ref_seq, $cds_cons_seq, $write_out);
                        my ($mut, $mutlines) = find_mutations($realn, 1, length($cds_string), $cds_protein, $ref_protein, $gene_name, $write_out);
                        report_mutations($mut, $mutlines, $cds_protein);

                        my $feat = new Bio::SeqFeature::Generic( -location => $splitlocation,
                                                                 -primary_tag => $feat_object->primary_tag,
                                                                 -tag => { 'gene' => $gene_name,
                                                                           'translation' => $cds_protein } );
                        $consensus->add_SeqFeature($feat);
                    }
                }
            } else {    # CDS not spliced
                $first_col = $aln->column_from_residue_number($ref, $feat_object->start);
                $last_col = $aln->column_from_residue_number($ref, $feat_object->end);

                $cds_aln = $aln->slice($first_col, $last_col, 0);
                if ($cds_aln->num_sequences != 2) {
                    say $outfile "no sequence";
                    next FEAT;
                }

                # extract consensus CDS and reference CDS sequences
                my $cds_seq = $cds_aln->get_seq_by_id($consensus_id);
                ($cds_string = $cds_seq->seq) =~ tr/\.//d;    # remove gaps from consensus


                #say $gene_name . "\t\t\t" . "coords: ", $feat_object->start, "-", $feat_object->end;
                my $cds_protein = translate_as_string($cds_string);
                #say 'cds nucleotides: ', $cds_string;
                #say 'ref nucleotides: ', $feat_object->spliced_seq->seq;
                #say 'cds protein: ', $cds_protein;
                #say 'ref protein: ', $ref_protein;

                # report mutations
                my $orig_start = rindex($clean_seq, $cds_string) + 1;
                my $orig_end = $orig_start + length($cds_string) - 1;
                my ($mut, $mutlines) = find_mutations($cds_aln, $orig_start, $orig_end, $cds_protein, $ref_protein, $gene_name, $write_out);
                report_mutations($mut, $mutlines, $cds_protein);

                # create new feature for CDS/gene and add to consensus sequence object
                my $feat = new Bio::SeqFeature::Generic( -start => $orig_start,
                                                         -end => $orig_end,
                                                         -primary_tag => $feat_object->primary_tag,
                                                         -tag => { 'gene' => $gene_name,
                                                                   'translation' => $cds_protein } );
                $consensus->add_SeqFeature($feat);
            }
            say $outfile "Length: " . length($cds_string) . " bp\t\tRef length: " . length($feat_object->spliced_seq->seq) . " bp";

            # add gene to concatednated seqequence genbank object
            if ($gene_name eq 'envelope') {
                $seqs{'env'} = $cds_string;
            }
            else {
                $seqs{$gene_name} = $cds_string;
            }


        }
    }
    # output annotated consensus sequence as genbank file
    $output = new Bio::SeqIO( -file => ">" . $outfilename . ".gb",
                              -format => "genbank" );
    $output->write_seq($consensus);
    p %seqs;

    #@gene_names = ('env');
    @gene_names = ('gag','pol','vif','vpr','tat','rev','vpu','env','nef');
    my $full_seq;
    my $prev_end = 1;
    for my $gene_name (@gene_names) {
        if (exists($seqs{$gene_name})){
            my $feat = new Bio::SeqFeature::Generic( -start => $prev_end,
                                                    -end => $prev_end + length($seqs{$gene_name}) - 1,
                                                    -primary_tag => 'CDS',
                                                    -tag => { 'gene' => $gene_name } );
            $concat->add_SeqFeature($feat);
            $prev_end = $prev_end + length($seqs{$gene_name});
            #say 'gene name' . $gene_name;
            #say 'seq' . $seqs{$gene_name};
            #say 'length seq' . length($seqs{$gene_name});
            #say 'prev end '.$prev_end;
            #$full_seq .= 'NNN' . $seqs{$gene_name};
            $full_seq .= $seqs{$gene_name};
            $concat->seq($full_seq);
        }
    }
    say $outseqfile $full_seq;

    # write out concatenated sequence of all regions
    $concat_output = new Bio::SeqIO( -file => ">" . $outfilename . "_concat.gb",
                              -format => "genbank" );
    $concat_output->write_seq($concat);

}


sub report_mutations {
    my ($mut, $mutlines, $cds_protein) = @_;
    my %mutations = %$mut;
    my @mutation_lines = @$mutlines;

    # determine if protein is productive or not
    my $productive = ( ($cds_protein =~ tr/\*//) == 1
                        && index($cds_protein, '*') >= ceil(length($cds_protein)*$min_length) )
                        ? 'yes' : 'no';
    say $outfile 'Productive: ' . $productive . "\t\tMismatches: " . $mutations{mismatches} .
                "\t\tInsertions: " . $mutations{insertions} . "\t\tDeletions: " . $mutations{deletions} . "\n";
    say $outfile 'Mutations:';

    my $tb = create_table();
    for (@mutation_lines) { $tb->load($_); }
    print $outfile $tb;
}


sub create_table {
    my $tb = Text::Table->new(
        { title => "\t\t\n\t\t",
           align => 'center',
           align_title => 'center',
        },
        { title => "Position\n--------",
           align => 'center',
           align_title => 'center',
        },
        { title => "RefBase\n-------",
           align => 'center',
           align_title => 'center',
        },
        { title => "Base\n----",
           align => 'center',
           align_title => 'center',
        },
        { title => "RefAA\n-----",
           align => 'center',
           align_title => 'center'
        },
        { title => "AA\n--",
           align => 'center',
           align_title => 'center'
        },
        { title => "Equiv\n-----",
           align => 'center',
           align_title => 'center'
        },
        { title => "FrameShift\n----------",
           align => 'center',
           align_title => 'center'
        },
    );
    return $tb;
}


main();
close $outfile;
