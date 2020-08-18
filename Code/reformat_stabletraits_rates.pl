#!/usr/bin/env perl 

# This Perl 5 script is meant to be called by the fit_stable_model.R
# script. It reformats the output of stabletraits, so that it is easier
# to analyse further.

use strict;
use warnings;

use feature qw(say);

use Statistics::R;

####################
# M A I N  C O D E #
####################

# Get the working directory.
my $working_dir = $ARGV[0];

# Figure out the number of species and their order in the phylogeny.
my $sp_number  = get_number_of_species($working_dir);
my $sp_indices = match_species_to_numbers($working_dir);

# Open the input file and the output file for writing.
open my $out_fh, '>', $working_dir . 'out.brlens_reformatted' or die $!;
open my $in_fh,  '<', $working_dir . 'out.brlens'             or die $!;

# For each line of the input file...
while ( my $line = <$in_fh> )
{

    # ... capture the node IDs and the rate of the branch.
    if ( $line =~ /^(\w+)->(\w+)\s.*?(\d+[.]\S+)$/ )
    {
        my $node_1 = $1;
        my $node_2 = $2;
        my $rate   = $3;

        if ( $node_1 =~ /^n(\d+)/ )
        {
            $node_1 = $1 + $sp_number + 1;
        }
        else
        {
            $node_1 = $sp_indices->{$node_1};
        }

        if ( $node_2 =~ /^n(\d+)/ )
        {
            $node_2 = $1 + $sp_number + 1;
        }

        # Print the node numbers and the rate.
        say {$out_fh} $node_2 . "\t" . $node_1 . "\t" . $rate;
    }
    else
    {

        # If we are at the very first line of the input file,
        # just print the header.
        say {$out_fh} "Node_1\tNode_2\tRate";
    }
}

close $in_fh;
close $out_fh;

###################################
# S  U  B  R  O  U  T  I  N  E  S #
###################################

# This subroutine counts the number of species in the current dataset.
sub get_number_of_species
{
    my ($dir) = @_;

    my $line_counter = 0;

    local $/ = "\n";
    open my $in_fh, '<', $working_dir . 'dataset.txt' or die $!;

    while ( my $line = <$in_fh> )
    {
        $line_counter++;
    }

    close $in_fh;

    return $line_counter;
}

# This function calls R and maps each species to its index number in
# the phylogeny.
sub match_species_to_numbers
{
    my ($dir) = @_;

    my $R = Statistics::R->new;

    $R->run(<<"END");
	library(ape)
	tree <- read.tree(paste("$dir", "tree.nwk", sep = ""))
END

    my $sp_list = $R->get('tree$tip.label');
    my %sp_indices;

    for ( 0 .. @{$sp_list} - 1 )
    {
        $sp_indices{ $sp_list->[$_] } = $_ + 1;
    }

    return \%sp_indices;
}
