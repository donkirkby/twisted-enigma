import typing
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from collections import defaultdict, Counter

# noinspection PyPackageRequirements
from Bio.Alphabet import IUPAC
# noinspection PyPackageRequirements
from Bio.Seq import Seq
# noinspection PyPackageRequirements
from Bio.Data.CodonTable import unambiguous_dna_by_name, CodonTable


def parse_args():
    parser = ArgumentParser(
        description='Find words that can be made from codons',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('dictionary',
                        type=FileType(),
                        nargs='?',
                        default='/usr/share/dict/british-english')
    return parser.parse_args()


def find_sequences(word: str,
                   reverse_word: str,
                   source_codons: typing.Mapping[str, set]) -> typing.List[str]:
    sequences = {''}
    for amino1, amino2 in zip(word, reversed(reverse_word)):
        new_sequences = set()
        for nucs1 in source_codons[amino1]:
            for nucs2 in source_codons[amino2]:
                rev_nucs2 = str(Seq(nucs2).reverse_complement())
                if nucs1 == rev_nucs2:
                    for old_seq in sequences:
                        new_sequences.add(old_seq + nucs1)
        sequences = new_sequences
    return list(sequences)


def filter_sequences_by_count(sequences: typing.List[str]) -> typing.List[str]:
    filtered_sequences = []
    for sequence in sequences:
        counts = Counter(sequence)
        if counts['G'] + counts['C'] > 10:
            continue
        if counts['A'] + counts['T'] > 10:
            continue
        filtered_sequences.append(sequence)
    return filtered_sequences


def main():
    args = parse_args()
    complementary_aminos = defaultdict(set)
    source_codons = defaultdict(set)
    table: CodonTable = unambiguous_dna_by_name['Standard']
    for codon, amino in table.forward_table.items():
        source_codons[amino].add(codon)
        rev_nucs = Seq(codon).reverse_complement()
        rev_amino = rev_nucs.translate()
        complementary_aminos[amino].add(rev_amino[0])

    for amino, rev_aminos in sorted(complementary_aminos.items()):
        print(f'{amino} <=> {", ".join(sorted(rev_aminos))}')

    clean_words = (word.upper().rstrip() for word in args.dictionary)
    sized_words = (word for word in clean_words if len(word) == 6)
    amino_letters = set(IUPAC.protein.letters)
    word_list = []
    for word in sized_words:
        extra_chars = set(word) - amino_letters
        if not extra_chars:
            word_list.append(word)
            if __name__ == '__live_coding__' and len(word_list) == 100:
                break
    all_words = set(word_list)
    word_list.sort()
    for word in word_list:
        reverse_words = {''}
        for amino in word:
            new_reverse_words = set()
            for reverse_amino in complementary_aminos[amino]:
                for reverse_word in reverse_words:
                    new_reverse_words.add(reverse_amino + reverse_word)
            reverse_words = new_reverse_words
        valid_words = reverse_words & all_words
        for reverse_word in valid_words:
            if reverse_word < word:
                continue
            seqs = find_sequences(word, reverse_word, source_codons)
            filtered_seqs = filter_sequences_by_count(seqs)
            if reverse_word == word and filtered_seqs:
                print('vvvvvvvvvvvvvvvvvv')
            for nucs in filtered_seqs:
                print(word, reverse_word, '<=>', nucs)

    seq = Seq('AGTGCTGATATCAGCACT')
    print(seq)
    print(seq.reverse_complement())
    print('Done.')


"""
ACIDIC ADIDAS <=> GCTTGCATCGATATCTGC
AFFIRM HANKER <=> GCGTTCTTTATTCGCATG
AFRICA CANTER <=> GCGTTCCGTATTTGCGCA
AGINGS TAINTS <=> GCTGGTATTAATGGCAGT
vvvvvvvvvvvvvvvvvv
AIDING AIDING <=> GCCATTGATATCAATGGC
AISLES GLEANS <=> GCTATTAGCCTCGAGTCC
ALASKA SLACKS <=> GCTCTTGCAAGCAAGGCT
ALASKA CLACKS <=> GCTCTTGCAAGCAAGGCA
ALIENS TILDES <=> GCTCTCATCGAGAATAGT
ALISSA CRANKS <=> GCTCTTATTAGCTCGGCA
ALISSA CRANKS <=> GCTCTTATTAGCTCTGCA
ALISSA CRANES <=> GCTCTCATTAGCTCTGCA
ALISSA CRANES <=> GCTCTCATTAGCTCGGCA
ALLAYS TICKER <=> GCGCTCCTTGCATATAGT
ALMACH MASHES <=> GCTCTCATGGCTTGCCAT
ALPINE FINGER <=> GCGCTCCCCATTAATGAA
ALPINE LINGER <=> GCGCTCCCCATTAATGAG
ALPINE FINGER <=> GCGCTCCCCATTAATGAA
ALPINE LINGER <=> GCGCTCCCCATTAATGAG
ALYSSA CARIES <=> GCTCTCTATTCTAGCGCA
ALYSSA CARIES <=> GCTCTCTATTCGAGCGCA
ALYSSA CRAVES <=> GCTCTCTACAGCTCTGCA
ALYSSA CARVES <=> GCTCTCTACTCTAGCGCA
ALYSSA STAVES <=> GCTCTCTACAGCAGTGCT
ANGELA CELTIC <=> GCAAATGGTGAGCTCGCA
ARLENE FILETS <=> GCTCGTCTCGAGAATGAA
ARLINE LINEAR <=> GCGCGCCTCATTAATGAG
ARMANI NIGHTS <=> GCTCGTATGGCCAATATT
ASIDES GLINTS <=> GCTAGTATTGATGAGTCC
ASLANT RISERS <=> GCTTCGCTCGCTAATACG
ASLANT RISERS <=> GCTTCTCTCGCTAATACG
ASSERT SPLATS <=> GCTAGTAGCGAGAGGACT
ASSIST RADARS <=> GCTTCTAGCATCAGCACG
CALIFS GENERA <=> TGCGCGCTCATTTTCTCC
CANALS RESIST <=> TGTGCTAATGCTCTCTCG
CANALS RESIST <=> TGTGCTAATGCTCTCTCT
CASALS RECAST <=> TGTGCTAGCGCACTCTCT
vvvvvvvvvvvvvvvvvv
CICADA CICADA <=> TGCATCTGCGCAGATGCA
CLAIRE LANCET <=> TGTCTCGCAATTCGCGAG
FAIRER PLANCK <=> TTTGCAATTCGCGAGAGG
FAIRER PLANCK <=> TTTGCAATTCGCGAGCGG
FERVID INHALE <=> TTCGAGCGCGTGATTGAT
FETISH MANGLE <=> TTCGAGACCATTAGCCAT
FIDDLE LEVINE <=> TTCATTGATGACCTCGAG
FINALS RESIDE <=> TTCATCAATGCTCTCTCG
GALENA RIFEST <=> GGTGCTCTCGAAAATGCG
GINGER PLAINT <=> GGTATTAATGGCGAGAGG
GINGER PLAINT <=> GGTATTAATGGCGAGCGG
GINGER PLAINS <=> GGAATTAATGGCGAGCGG
GINGER PLAINS <=> GGAATTAATGGCGAGAGG
GINGER PLAIDS <=> GGAATCAATGGCGAGAGG
GINGER PLAINT <=> GGTATTAATGGCGAGAGG
GINGER PLAINT <=> GGTATTAATGGCGAGCGG
GINGER PLAINS <=> GGAATTAATGGCGAGCGG
GINGER PLAINS <=> GGAATTAATGGCGAGAGG
GINGER PLAIDS <=> GGAATCAATGGCGAGAGG
GIRDER PLIANT <=> GGTATTCGCGATGAGAGG
GLASER SLACKS <=> GGACTTGCAAGCGAGAGA
GLIDER SLIDES <=> GGACTCATCGATGAGAGA
GLIDER SLIDES <=> GGACTCATCGATGAGCGA
GLIDER SLINKS <=> GGACTTATTGATGAGCGA
vvvvvvvvvvvvvvvvvv
GLIDES GLIDES <=> GGACTCATCGATGAGTCC
GLITCH MAGNET <=> GGTCTCATTACCTGCCAT
ISLAMS THREAD <=> ATCAGCCTCGCGATGAGT
ITASCA STARRY <=> ATAACGGCGAGCTGTGCT
LANCER PLAICE <=> CTCGCAAATTGCGAGAGG
MANICS RADISH <=> ATGGCTAATATCTGCTCG
MARCIA SNATCH <=> ATGGCACGTTGCATTGCT
PLAIDS RINGER <=> CCTCTCGCCATTGATTCT
PLAIDS RINGER <=> CCGCTCGCCATTGATTCT
PLAIDS RINGER <=> CCTCTCGCCATTGATTCG
PLAINS RINGER <=> CCTCTCGCCATTAATTCT
PLAINS RINGER <=> CCGCTCGCCATTAATTCT
PLAINS RINGER <=> CCTCTCGCCATTAATTCG
PLAINS RINGER <=> CCGCTCGCCATTAATTCG
PLAINT SINGER <=> CCGCTCGCCATTAATACT
PLAINT SINGER <=> CCTCTCGCCATTAATACT
PLAINT RINGER <=> CCGCTCGCCATTAATACG
PLAINT RINGER <=> CCTCTCGCCATTAATACG
PLANET SLICER <=> CCGCTCGCAAATGAGACT
PLANET SLICER <=> CCTCTCGCAAATGAGACT
PLAYER SLICER <=> CCTCTCGCATATGAGCGA
PLAYER SLICER <=> CCGCTCGCATATGAGAGA
PLAYER SLICER <=> CCTCTCGCATATGAGAGA
PLIANT SIGNER <=> CCGCTCATTGCCAATACT
PLIANT SIGNER <=> CCTCTCATTGCCAATACT
PLYING TIDIER <=> CCGCTCTATATCAATGGT
PRAWNS RIPSAW <=> CCACGCGCTTGGAATTCT
RAFFIA SNEERS <=> AGAGCGTTCTTCATTGCT
RAFFIA SNEERS <=> CGAGCGTTCTTCATTGCT
RAIDER SLINGS <=> AGAGCCATTGATGAGAGA
RAIDER SLINGS <=> CGAGCCATTGATGAGCGA
RAIDER SLINGS <=> CGAGCCATTGATGAGAGA
RAIDER SLINGS <=> AGAGCCATTGATGAGCGA
vvvvvvvvvvvvvvvvvv
SADIST SADIST <=> AGTGCTGATATCAGCACT
vvvvvvvvvvvvvvvvvv
SAFEST SAFEST <=> AGTGCTTTCGAAAGCACT
SALIVA SHYEST <=> AGTGCTCTCATAGTGGCT
SALVER SLYEST <=> AGTGCTCTCGTAGAGCGA
SALVER SLYEST <=> AGTGCTCTCGTAGAGAGA
SLAYER SLICER <=> TCGCTCGCATATGAGAGA
SLAYER SLICER <=> TCTCTCGCATATGAGCGA
SLAYER SLICER <=> TCTCTCGCATATGAGAGA
SLAYER SLICER <=> TCGCTCGCATATGAGCGA
vvvvvvvvvvvvvvvvvv
SLIDER SLIDER <=> TCTCTCATCGATGAGAGA
SLIDER SLIDER <=> TCGCTCATCGATGAGAGA
SLIDER SLIDER <=> TCTCTCATCGATGAGCGA
SLIDER SLIDER <=> TCGCTCATCGATGAGCGA
vvvvvvvvvvvvvvvvvv
TRACTS TRACTS <=> ACTCGTGCATGCACGAGT
TRAITS TRYSTS <=> ACTCGTGCTATAACGAGT
vvvvvvvvvvvvvvvvvv
VANISH VANISH <=> GTGGCTAATATTAGCCAC
VASTLY VESTRY <=> GTAGCGAGTACTCTCTAC
"""

main()
