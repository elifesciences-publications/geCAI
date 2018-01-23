#! /usr/local/bin/perl -w

use Statistics::RankCorrelation;
use Getopt::Std;

sub set_starting_scores	{
						%START = ();
						%synon = (TTT, F, TCT, S, TAT, Y, TGT, C, TTC, F, TCC, S, TAC, Y, TGC, C, TTA, L, TCA, S, TTG, L, TCG, S, TGG, W, CTT, L, CCT, P, CAT, H, CGT, R, CTC, L, CCC, P, CAC, H, CGC, R, CTA, L, CCA, P, CAA, Q, CGA, R, CTG, L, CCG, P, CAG, Q, CGG, R, ATT, I, ACT, T, AAT, N, AGT, S, ATC, I, ACC, T, AAC, N, AGC, S, ATA, I, ACA, T, AAA, K, AGA, R, ATG, M, ACG, T, AAG, K, AGG, R, GTT, V, GCT, A, GAT, D, GGT, G, GTC, V, GCC, A, GAC, D, GGC, G, GTA, V, GCA, A, GAA, E, GGA, G, GTG, V, GCG, A, GAG, E, GGG, G);
						%aa_groups = ();						
						foreach my $c (keys(%synon))
							{
							my $a = $synon{$c};
							$aa_groups{$a}{$c} = 1;
							}
						
						if (($weights_start == 1)&&(-e $weights_file))
							{
							open FILE, $weights_file;
							while (<FILE>)
								{
								my $l = $_;
								chomp $l;
								my @t = split(/\t/, $l);
								$START{$t[0]} = $t[1];
								}
							close FILE;
							}						
						else
							{
							foreach my $e (keys(%synon))
								{
								$START{$e} = rand();
								}							
							}
						
						my @sort = sort {$START{$a} <=> $START{$b}} keys %START;
						foreach my $e (@sort)
							{								
							$START{$e} = $START{$e}/$START{$sort[-1]};								
							}
						}
						
sub read_sequences		{						
						my $file = $_[0];
						chomp $file;
						open FILE, $file;
						my @in = <FILE>;
						close FILE;
						%seqs = ();
						my $j = join '', @in;
						$j =~ s/>/\n\n>/gmis;
						$j = $j."\n\n\n";						
						while ($j =~ m/>(.+?)\n(.+?)\n\n/gmis)
							{
							my %local = ();
							my $a = $1;							
							my $s = $2;							
							$s =~ tr/a-z/A-Z/;
							$s =~ s/[\s\n]//gmis;
							if ((exists $exp{$a})&&(length $s > 300))
								{								
								unless (($s =~ m/N/)||($s =~ m/Y/)||($s =~ m/R/)||($s =~ m/K/))
									{
									$seqs{$a} = $s;									
									push(@order, $a);																	
									while ($s =~ m/(\w{3})/gmis)
										{
										my $w1 = $1;
										unless (($w1 eq 'TGA')||($w1 eq 'TAG')||($w1 eq 'TAA'))
											{
											++$codons{$a}{$w1};
											++$codon_count{$a};
											++$count{$w1};
											if (exists $choice{$w1})
												{
												++$local{$w1};
												}
											}
										}									
									}
								}							
							}					
						}
						
sub mutate	{
			my %temp = %BEST;
			my %local = ();
			my @tmp = keys %temp;
			my $mut = int(rand($#tmp - 1));
			my $c = 0;
			while ($c < $mut)
				{				
				my $p = int(rand($#tmp + 1));			
				unless (exists $local{$p})
					{
					my $v = $temp{$tmp[$p]};	
					my $x = rand(0.1) + $min;
					my $s = rand(1);
					if ($s > 0.5)
						{
						$temp{$tmp[$p]} += $x;
						if ($temp{$tmp[$p]} > 1)
							{
							$temp{$tmp[$p]} = 1;
							}					
						}
					else
						{
						$temp{$tmp[$p]} -= $x;
						if ($temp{$tmp[$p]} < $min)
							{
							$temp{$tmp[$p]} = $min;
							}					
						}
					unless ($v == $temp{$tmp[$p]})
						{
						$local{$p} = 1;
						++$c;
						}
					}
				}			
			my @sort = sort {$temp{$a} <=> $temp{$b}} keys %temp;
			$temp{$sort[-1]} = 1;
			%CURRENT = %temp;
			foreach my $a (keys(%CURRENT))
				{
				$LOG{$a} = log($CURRENT{$a});
				}				
			}
			
sub mutate1	{			
			my %temp = %BEST;
			my @tmp = keys %temp;
			my $pass = 0;
			while ($pass == 0)
				{
				my $p = int(rand($#tmp + 1));
				my $v = $temp{$tmp[$p]};			
				my $x = rand(0.5);
				my $s = rand(1);
				if ($s > 0.5)
					{
					$temp{$tmp[$p]} += $x;
					if ($temp{$tmp[$p]} > 1)
						{
						$temp{$tmp[$p]} = 1;
						}					
					}
				else
					{
					$temp{$tmp[$p]} -= $x;
					if ($temp{$tmp[$p]} < $min)
						{
						$temp{$tmp[$p]} = $min;
						}					
					}
				unless ($v == $temp{$tmp[$p]})
					{
					$pass = 1;
					}					
				}		
			
			my @sort = sort {$temp{$a} <=> $temp{$b}} keys %temp;
			$temp{$sort[-1]} = 1;
			%CURRENT = %temp;	
			foreach my $a (keys(%CURRENT))
				{
				$LOG{$a} = log($CURRENT{$a});
				}				
			}
			
sub compute_CAI_fast	{
						my $a = $_[0];												
						my $value = 0;						
						foreach my $c (keys(%{$codons{$a}}))
							{										
							$value += ($codons{$a}{$c})*($LOG{$c});								
							}
						$value = exp($value/$codon_count{$a});											
						return $value;				
						}
				
sub	score_all_genes	{									
					%SCORE = ();
					@set1 = ();					
					foreach my $a (@order)
						{						
						$SCORE{$a} = &compute_CAI_fast($a);							
						push (@set1, $SCORE{$a});
						}
					}
					
sub assess_update	{
					my $current = $_[0];
					my $previous = $winning_round;					
					
					my $c1 = Statistics::RankCorrelation->new(\@set1, \@set2);
					my $current_cor = $c1->spearman;					
																																					
					if ($current_cor > $previous_cor)
						{
						%BEST = %CURRENT;						
						$winning_round = $current;
						$previous_cor = $current_cor;						
						open OUT, ">geCAI_codon_weights.txt";
						foreach my $a (keys(%BEST))
							{
							print OUT "$a\t$BEST{$a}\n";
							}
						close OUT;
						}					
					my $vtem = $current/$print_freq;
					if ((int($vtem) == $vtem)||($current == 1))
						{
						print " $current\t$previous_cor\n";
						}							
					}
					
					
sub pearsons	{							
				my $i = 0;
				while ($i <= $#set1)
					{
					$sum_x += $set1[$i];
					$sum_y += $set2[$i];
					$sum_xy += ($set1[$i])*($set2[$i]);
					$sum_x2 += $set1[$i]**2;
					$sum_y2 += $set2[$i]**2;						
					++$i;
					}
																
				my $count = $#set1 + 1;
				my $var = ($sum_x2 - (($sum_x**2)/$i))/($i - 1);						
				my $cor = ((($count)*($sum_xy))-($sum_x)*($sum_y))/((sqrt(($count*$sum_x2 - ($sum_x**2))))*(sqrt(($count*$sum_y2 - ($sum_y**2)))));
				return $cor, $var;	
				}
					
sub def_param	{
				%synon = (TTT, F, TCT, S, TAT, Y, TGT, C, TTC, F, TCC, S, TAC, Y, TGC, C, TTA, L, TCA, S, TTG, L, TCG, S, TGG, W, CTT, L, CCT, P, CAT, H, CGT, R, CTC, L, CCC, P, CAC, H, CGC, R, CTA, L, CCA, P, CAA, Q, CGA, R, CTG, L, CCG, P, CAG, Q, CGG, R, ATT, I, ACT, T, AAT, N, AGT, S, ATC, I, ACC, T, AAC, N, AGC, S, ATA, I, ACA, T, AAA, K, AGA, R, ATG, M, ACG, T, AAG, K, AGG, R, GTT, V, GCT, A, GAT, D, GGT, G, GTC, V, GCC, A, GAC, D, GGC, G, GTA, V, GCA, A, GAA, E, GGA, G, GTG, V, GCG, A, GAG, E, GGG, G);
				foreach my $k (keys(%synon))
					{
					$count{$k} = 0;
					}
				@trans = (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y);
				foreach my $entry (@trans)
					{
					@$entry = ();
					}
				foreach my $codon (keys(%synon))
					{
					my $temp = $synon{$codon};
					push (@$temp, $codon);
					}
			
				%choice = qw/ATT 0 CTT 0 GTT 0 TGT 0 CCT 0 ACT 0 TCT 0 CGT 0/;
				@choice = sort {$a cmp $b} keys %choice;
				}					
							
sub read_expression_data		{
								%exp = ();
								open FILE, $exp_file;
								while (<FILE>)
									{
									my $l = $_;
									chomp $l;
									my @t = split(/\t/, $l);
									if ((exists $t[1])&&($t[1] > 1))
										{
										$exp{$t[0]} = $t[1];								
										}
									}
								close FILE;							
								}
						
sub run_program	{
				%codons = ();
				&def_param;
				&get_options;
				print "Generating random codon weights...\n";
				&set_starting_scores;	
				print "Reading expression data...\n";			
				&read_expression_data;				
				print "Reading sequence data...\n";				
				&read_sequences($seq_file);									
				print "Running MCMC...\n\n";
				%BEST = %START;
				%CURRENT = %START;
				foreach my $a (keys(%CURRENT))
					{
					$LOG{$a} = log($CURRENT{$a});
					}
				
				%SCORE = ();
				@set1 = ();
				@set2 = ();		
				
				&score_all_genes(0);							
				$winning_round = 0;			
						
				foreach my $a (@order)
					{					
					push (@set2, $exp{$a});								
					}
								
				my $c1 = Statistics::RankCorrelation->new(\@set1, \@set2);
				my $spearman = $c1->spearman;							
				my ($pearson, $var) = &pearsons(@set1);	
				$previous_cor = $spearman;								
				
				print " Gen.\tSpearman correlation\n";						
				my $i = 1;
				$tc = 0;
				while ($i <= $ngen)
					{
					if ($tc == 1)
						{
						&mutate1(\%BEST);
						$tc = 0;
						}
					else
						{
						&mutate(\%BEST);
						$tc = 1;
						}						
					&score_all_genes($i);
					&assess_update($i);
					++$i;
					}
				print "\nRUN COMPLETE:\n";
				print " Resulting geCAI codon weights can be found in geCAI_codon_weights.txt\n";
				&end_of_run;				
				}
				
sub exit_screen	{
				system "clear";
				print "\ngeCAI calculator version 1.0 Copyright (c) 2017 Steven Kelly\n";
				print "\nUSAGE:\n";
				print " geCAI.pl -i <SEQUENCE_FILE> -e <EXPRESSION_FILE>\n";
				print "\nOPTIONS:\n";
				print " -i <file>	FASTA format file containing multiple CDS sequences\n";
				print " -e <file>	Tab-delimitted text file of abundance estimates for CDS\n";
				print " -g <int>	Number of generations for MCMC [Default = 5000]\n";
				print " -p <int>	Print frequency [Default = 100]\n";
				print " -w <file>	Specify staring codon weights [Defualt = random]\n";
				&end_of_run;
				}
				
sub end_of_run	{				
				print "\nLICENSE:\n";
				print " Distributed under the GNU General Public License (GPLv3). See License.md\n";
				print "\nCITATION:\n";
				print " When publishing work that uses geCAI please cite\:\n";
				print " Nascimento & Kelly et al. (2018)\n";
				print "\n\n";
				exit;
				}

sub get_options	{
				&getopts('i:w:e:g:p:',\%parameters);
				if (exists $parameters{'i'})
					{
					$seq_file = $parameters{'i'};
					chomp $seq_file;					
					}				
				else
					{
					&exit_screen;					
					}
				if (exists $parameters{'e'})
					{
					$exp_file = $parameters{'e'};
					chomp $exp_file;					
					}				
				else
					{
					&exit_screen;					
					}
				if (exists $parameters{'w'})
					{
					$weights_file = $parameters{'w'};
					chomp $weights_file;
					$weights_start = 1;					
					}				
				else
					{
					$weights_start = 0;
					}					
				if (exists $parameters{'g'})
					{
					$ngen = $parameters{'g'};
					chomp $ngen;								
					}				
				else
					{
					$ngen = 5000;
					}
				if (exists $parameters{'p'})
					{
					$print_freq = $parameters{'p'};
					chomp $print_freq;								
					}				
				else
					{
					$print_freq = 100;
					}
				system "clear";
				print "\ngeCAI calculator version 1.0 Copyright (c) 2017 Steven Kelly\n";
				print "\n";
				
				$min = 0.01;						
				}
				
&run_program;

exit;


