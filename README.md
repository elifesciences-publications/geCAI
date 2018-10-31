

### This code is associated with the paper from de Freitas Nascimento et al., "Codon choice directs constitutive mRNA levels in trypanosomes". eLife, 2018. http://dx.doi.org/10.7554/eLife.32467


# geCAI
MCMC method for estimating geCAI codon weights from CDS and mRNA abundance data 

## Running geCAI
```
geCAI.pl -i <SEQUENCE_FILE> -e <EXPRESSION_FILE>

OPTIONS:
 -i <file>	FASTA format file containing multiple CDS sequences
 -e <file>	Tab-delimitted text file of abundance estimates for CDS
 -g <int>	Number of generations for MCMC [Default = 5000]
 -p <int>	Print frequency [Default = 100]
 -w <file>	Specify staring codon weights [Defualt = random]
```

## EXAMPLE FILES
Package contains example sequence and mRNA files. To run geCAI on these files execute
```
geCAI.pl -i Example_Sequences.fasta -e Example_Expression_Data.txt
```
## LICENSE:
 Distributed under the GNU General Public License (GPLv3). See License.md

## CITATION:
 When publishing work that uses geCAI please cite:
 Nascimento & Kelly et al. (2018)
