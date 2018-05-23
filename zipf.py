def rep(x, y):  # too lazy to go through and change R code + regex was being wonky -> 'rep\("([A-Z])", {1}([0-9])\)'
    return list(x * y)


codons = ["GCU", "GCC", "GCA", "GCG", "CGU", "CGC", "CGA", "CGG", "AGA", "AGG", "AAU", "AAC", "GAU", "GAC", "UGU",
          "UGC", "CAA", "CAG", "GAA", "GAG", "GGU", "GGC", "GGA", "GGG", "CAU", "CAC", "AUU", "AUC", "AUA", "UUA",
          "UUG", "CUU", "CUC", "CUA", "CUG", "AAA", "AAG", "AUG", "UUU", "UUC", "CCU", "CCC", "CCA", "CCG", "UCU",
          "UCC", "UCA", "UCG", "AGU", "AGC", "ACU", "ACC", "ACA", "ACG", "UGG", "UAU", "UAC", "GUU", "GUC", "GUA",
          "GUG", "UAA", "UGA", "UAG"]  # 64 codons
amino_acids = [item for sublist in  # flatten this boi into a list of 64 amino acids 1-letter abbreviations
               [rep("A", 4), rep("R", 6), rep("N", 2), rep("D", 2), rep("C", 2), rep("Q", 2), rep("E", 2), rep("G", 4),
                rep("H", 2), rep("I", 3), rep("L", 6), rep("K", 2), rep("M", 1), rep("F", 2), rep("P", 4), rep("S", 6),
                rep("T", 4), rep("W", 1), rep("Y", 2), rep("V", 4), rep("*", 3)] for item in sublist]
codon_table = dict((k, v) for k, v in zip(codons, amino_acids))


def translate(nucleotides):
    """
    Translates DNA to RNA
    :param nucleotides: the DNA sequence
    :return: an RNA sequence
    """
    nucleotides = nucleotides.replace("C", "X")
    nucleotides = nucleotides.replace("G", "C")
    nucleotides = nucleotides.replace("A", "U")
    nucleotides = nucleotides.replace("T", "A")
    nucleotides = nucleotides.replace("X", "G")
    return nucleotides.strip()  # newline removal


def ribosome(nucleotides, simplify):
    """
    Attempts to translate a nucleotide sequence, but
    :param nucleotides: a DNA nucleotide sequence
    :param simplify: if true, when a sequence isn't a multiple of 3, it will condense codons into a proetin
    :return: either a RNA sequence or a protein sequence
    """
    if len(nucleotides) % 3 == 0:
        return "".join([codon_table[nucleotides[i:i + 3]] for i in range(0, len(nucleotides), 3)])
    else:
        if simplify:
            usable_length = len(nucleotides) - (len(nucleotides) % 3)
            print nucleotides
            return "".join([codon_table[nucleotides[i:i + 3]] for i in range(0, usable_length, 3)]) + \
                   nucleotides[-(len(nucleotides) % 3):]
        else:
            return nucleotides


def unit_generator(input_file, size, simplify):
    """
    Returns a semantic unit of a given size. Multiples of 3 will be protein sequences
    :param input_file: the location of the file that contains the sequence to be digested
    :param size: the size in nucleotides of the desired genetic unit
    :param simplify: if true, when a sequence isn't a multiple of 3, it will condense codons into a protein
    :return: a generator to supply genetic units from the sequence
    """
    if size < 1 or int(size) != size:
        raise ValueError("Size must be positive integer")
    with open(input_file) as f:
        curr = f.read(size)
        while curr is not None:  # gradually read whole sequence in chunks of a specified size
            seq = translate(curr)
            yield ribosome(seq, simplify)
            curr = f.read(size)


def make_frequency_table(generator):
    """
    Creates a table of the frequencies of certain sized chunks in a genome
    :param input_file: the location of the file that contains the sequence to be digested
    :param size: the size in nucleotides of the desired genetic unit
    :param simplify: if true, when a sequence isn't a multiple of 3, it will condense codons into a protein
    :return: a dictionary containing the table of frequencies
    """
    gen = generator  # default for straight DNA sequence
    freq_table = dict()
    unit = next(gen)
    counter = 0
    while unit is not "":  # form table from units
        if unit not in freq_table:
            freq_table[unit] = 0
        freq_table[unit] += 1
        unit = next(gen)
        counter += 1
    s_t = sorted([(k, v) for k, v in freq_table.items()], key=lambda x: x[1], reverse=True)  # sort by freq, descending
    return s_t, counter


###########################################################################################
## DNA SEQUENCE
# frequencies = [make_frequency_table(unit_generator(input_file="nucs.fasta", size=siz, simplify=True)) for siz in sizes]

# PROTEIN SEQUENCE
input_file = "/Users/matthew/protein.fa"
with open(input_file) as f:
    body = f.readlines()
    amino_acids = []
    for line in body:
        if ">" in line:
            continue
        amino_acids.append(line.strip())
    amino_acids = "".join(amino_acids)

size = 4  # 4 amino acids
freq_table = dict()
counter = 0
for s in range(size):
    for i in range(0, len(amino_acids), size):  # form table from units
        unit = amino_acids[i:i + size]
        if unit not in freq_table:
            freq_table[unit] = 0
        freq_table[unit] += 1
        counter += 1
frequencies = [sorted([(k, v) for k, v in freq_table.items()], key=lambda x: x[1], reverse=True)]  # sort by descending
with open("zipf.freq.txt", "wb") as o: # write frequencies to a file
    for y in frequencies[0]:
        o.write("{}\n".format(y[1]))

### RNA SEQUENCE
# input_file = "rna.fa"
# with open(input_file) as f:
#     body = f.readlines()
# dna = []
# for line in body:
#     if ">" in line:
#         continue
#     dna.append(line.strip())
# amino_acids = ribosome(translate("".join(dna)), simplify=True)
# print amino_acids[:20]
# size = 4
# freq_table = dict()
# counter = 0
# for s in range(size):
#     for i in range(0, len(amino_acids), size):  # form table from units
#         unit = amino_acids[i:i + size]
#         if unit not in freq_table:
#             freq_table[unit] = 0
#         freq_table[unit] += 1
#         counter += 1
# frequencies = [sorted([(k, v) for k, v in freq_table.items()], key=lambda x: x[1], reverse=True)]  # sort by descending

for table in frequencies:  # write zipf enumeration to a file
    rel_freq = list(map(lambda t: float(t[1]) / float(counter), table))
    top_entry = rel_freq[0]
    zipf = list(map(lambda f: (f / top_entry) ** -1, rel_freq))

    # pritn first 100 to 200 entries for each step in teh calculation
    print("table:", table[:100])
    print("relative frequency:", rel_freq[:100])
    print("zipf:", zipf[:200])
    with open("zipf", "wb") as o:
        for y in zipf:
            o.write("{}\n".format(y))

### R code
## for random sampling:
#     writeLines(paste(sample(nuc, 120000, replace = T), collapse=""), "/Users/matthew/PycharmProjects/BigBoiCode/nucs.fasta")
## for plotting:
#     library(tidyverse)
#     par(mfrow=c(2,2))
#     data <- read_csv("/Users/matthew/PycharmProjects/BigBoiCode/zipf", col_names = F)
#     `Frequency Ranking` <-  (1:length(data[[1]]))
#     `Zipf Enumeration` <- (data[[1]])
#     plot(`Frequency Ranking`, `Zipf Enumeration`, type="l")
#     title("Evaluating for Zipf Distribution")

#     `log(Frequency Ranking)` <-  log(1:length(data[[1]]))
#     `log(Zipf Enumeration)` <- log(data[[1]])
#     plot(`log(Frequency Ranking)`, `log(Zipf Enumeration)`, type="l")
#     title("Evaluating for Zipf Distribution (logarithmic)")

#     data <- read_csv("/Users/matthew/PycharmProjects/BigBoiCode/zipf.freq.txt", col_names = F)
#     `Rank` <-  (1:length(data[[1]]))
#     `Frequency` <- (data[[1]])
#     plot(`Rank`, `Frequency`, type="l")
#     title("Evaluating Frequency vs Rank");data <- read_csv("/Users/matthew/PycharmProjects/BigBoiCode/zipf.freq.txt", col_names = F)

#     `log(Rank)` <-  log(1:length(data[[1]]))
#     `log(Frequency)` <- log(data[[1]])
#     plot(`log(Rank)`, `log(Frequency)`, type="l")
#     title("Evaluating Frequency vs Rank (logarithmic)")
