# Twisted Enigma
This explains how I created the [Twisted Enigma] puzzle that I posted on Stack
Exchange. Obviously, it's a complete spoiler for the puzzle, so go solve it
first, if you don't want it spoiled.

Still here? Then here's the puzzle: what is special about this?

![puzzle photo]

I gave some extra clues, but the main objective was to identify what was in the
photo, and what was special about it. I was looking for answers at five levels:

1. It's DNA. In the first day after I published, four out of five answers
    identified it as DNA.
2. There are six codons in the sequence. If you've forgotten your high school
    biology, DNA is made up of two chains of nucleotides, and codons are groups
    of three nucleotides. Each codon is a code for an amino acid, and each
    amino acid is represented by a letter. The sequence in the picture could
    translate to a six-letter word. One answer and one comment got this far.
    (Please don't put spoilers in comments, but well done anyway.)
3. The colour scheme is the closest I could find to a standard: the [DRuMS]
    colour scheme. The nucleotides A, C, T, G are represented by blue, red,
    yellow, and green with the mnemonics Azure, Crimson, Tweety bird, and Green.
    Together with the [codon table] you could then decode the sequence that
    starts with yellow, red, green as "SLIDER". That matches a couple of the
    clues to let solvers know they're on the right track. One answer got this
    far, although it didn't use the same colour scheme. It used the
    complementary colour scheme, so it just read the other chain. I tried to
    frame the photo so that the chain I wanted solvers to tackle first was
    unobstructed and read top to bottom, but it was impossible to be clear
    without some way to identify the 5' and 3' ends of the molecule.
4. Now it gets more subjective. What I was looking for was for the solver to
    notice that the two chains encoded "SLIDER" in opposite directions. One
    answer got this.
5. The final level was to notice that although the two chains encoded "SLIDER",
    the actual nucleotides were different. Each amino acid can be encoded by
    more than one set of nucleotides. One answer got this complete solution.
    Well done, [tmpearce].

If you're interested in how I constructed the puzzle, I had the idea that I
could find an interesting genetic sequence and build it with a game that I
recently found at a thrift store: Splice. I decided that encoding a word would
be fun, and then the idea of palindromes restricted the search space.

I wrote the [`puzzle.py`] script to search through the dictionary, and look for
sequences that encode a word in both directions. It highlights sequences where
both directions encode the same word. After I did that and found that there
were several words that worked, I had the idea to look for sequences where the
amino acid sequences were the same in both directions, but the nucleotide
sequences weren't. SLIDER was the only word in my dictionary that worked, so
I decided to go with that.

Once I had the sequence I wanted, I discovered a couple of problems with my
Splice game. First, the helix is left handed, while most DNA molecules are right
handed. Luckily, I could fix that easily by flipping the photograph. Harder to
fix was the fact that the nucleotide pairs in the game didn't connect the
colours according to the standard colour scheme. I think it's very strange that
they used the standard colours, but didn't pair them correctly. Considering that
and the wrong handedness, I would guess that they didn't ask any biologists
before they manufactured the game. I briefly considered taking apart the pairs
and regluing them, but instead I decided to fix it by adjusting the colours in
the photo. With some good advice from Mike O'Shaughnessy, I wrote the
[`swap_colours.py`] script to swap the blue and red, as well as flip the photo.
Hopefully, nobody noticed the photo manipulation, but if you look at the edges
of the blue and red sections, you can see it. Here's the original photo:

![Original photo]

That's all the gory details. Send me a tweet @donkirkby or open an issue if you
can find other words that work, or if you'd like to collaborate on designing a
puzzle. If you enjoyed this, you might like my [Donimoes] collection.

[Twisted Enigma]: https://puzzling.stackexchange.com/q/91664/38
[puzzle photo]: swapped.jpg
[DRuMS]: https://proteopedia.org/wiki/index.php/DRuMS#Additional_DRuMS_ColorKey_Templates
[codon table]: https://en.wikipedia.org/wiki/DNA_codon_table
[tmpearce]: https://puzzling.stackexchange.com/a/91687/38
[`puzzle.py`]: puzzle.py
[`swap_colours.py`]: swap_colours.py
[Original photo]: helix.jpg
[Donimoes]: https://donkirkby.github.io/donimoes/
