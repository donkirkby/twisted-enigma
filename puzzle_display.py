# noinspection PyPackageRequirements
from Bio.Seq import Seq


def main():
    puzzle = 'GTGGCTAATATTAGCCAC'
    seq = Seq(puzzle)
    aminos = seq.translate()
    rev_seq = seq.reverse_complement()
    rev_aminos = rev_seq.translate()
    display_colours = {'A': 'Blue',  # Azure
                       'C': 'Red',  # Crimson
                       'G': 'Green',
                       'T': 'Yellow'}  # Tweety Bird
    # Swap red and blue, because the game is wrong, and we'll adjust the photo.
    display_colours['A'], display_colours['C'] = \
        display_colours['C'], display_colours['A']
    for nuc in seq:
        print(display_colours[nuc])
    print(f'Sequence: {seq}')
    print(f'Aminos: {aminos}')
    print(f'Reverse Sequence: {rev_seq}')
    print(f'Reverse Aminos: {rev_aminos}')


main()
