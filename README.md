# Genomics & Natural Language Processing
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The purpose of this project is to determine if DNA/RNA or protein sequences hold some of the same properties as text. Much like language has grammar and syntax, genomic data may have rules for combining meaningful units of information. Sequences outline for the cell how to assemble and fold a chain of amino acids into a protein, just as strings of words woven together create a book to be read — the pages of the book have to be meaningful as well, as a malformed protein will be useless to the cell.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This analagous nature suggests there may be value in analyzing genomic information using Natural Language Processing techniques. As a personal experiment to further prove the validity of this analogy, I decided to investigate how well genomic information could be conformed to a Zipfian distribution as a sort of litmus test for compatibility. Literature ranging from _Beowulf_ to hieroglypics to mail order catalogs will show up with evidence of Zipf's Law. It is bounded neither by language nor media, with the former especially evidenced in an analysis of Wikipedia text:
![Rank versus frequency for the first 10 million words in 30 Wikipedias (dumps from October 2015) in a log-log scale](/images/zipf_wikipedia.png?raw=true)
## Zipf's Law
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Zipf's law states that "given some corpus of natural language utterances, the frequency of any word is inversely proportional to its rank in the frequency table. Thus the most frequent word will occur approximately twice as often as the second most frequent word, three times as often as the third most frequent word." Mathematically,  
* _N_ as the number of elements
* _k_ as their rank
* _s_ as the value of the exponent characterizing the distribution  
&nbsp;
![Zipf Equation](/images/equation.png?raw=true)
&nbsp;  
## Results
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This in mind, I sourced a sequence from [_Bos Taurus_](ftp://ftp.ncbi.nih.gov/genomes/Bos_taurus/protein/protein.fa.gz) to analyze, with the semantic unit being an amino acid. This _corpus_ would end up being ~ 50 MB. I did write code to work with DNA and RNA as well, with the only difference being DNA would have to be transcribed into RNA and then translated into amino acids. The redundancy of the genetic code means that the same protein can have multiple different DNA/RNA translations (up to 6, in some cases), but this was disregarded as the cell only cares for what the amino acid in the end, completely ignoring any silent mutations.
![](/images/zipf_and_freq.png?raw=True)
In the upper right figure, it can be seen that the line nearly straightens out, but does not become linear. Ideally, the relationship between the Zipf Enumeration (given by the formula above) and the ranking of the word by frequency would be 1:1, but this is clearly not the case. Instead, it would appear that there are a number of increasingly rare "words", throwing off this transformation. However, extending this code to more and more protein sequences, I suspect the average slope would approach 1. 
