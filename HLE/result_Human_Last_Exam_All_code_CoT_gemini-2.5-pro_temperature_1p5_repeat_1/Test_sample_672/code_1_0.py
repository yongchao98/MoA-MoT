import sys

def perform_scansion():
    """
    This function determines the scansion of the given line of poetry
    and prints the result using 'x' for unstressed and '/' for stressed syllables.

    The line is: "And in the letter, my cousin mentions a piece of advice"

    Syllable breakdown and stress analysis:
    1.  And: Unstressed (x)
    2.  in: Unstressed (x)
    3.  the: Unstressed (x)
    4.  let-ter: Stressed, Unstressed (/x)
    5.  my: Unstressed (x)
    6.  cou-sin: Stressed, Unstressed (/x)
    7.  men-tions: Stressed, Unstressed (/x)
    8.  a: Unstressed (x)
    9.  piece: Stressed (/)
    10. of: Unstressed (x)
    11. ad-vice: Unstressed, Stressed (x/)
    
    The resulting pattern is a concatenation of these markers.
    """
    
    # Concatenating the stress markers for each syllable
    scansion_result = "x" + "x" + "x" + "/" + "x" + "x" + "/" + "x" + "/" + "x" + "x" + "/" + "x" + "x" + "/"
    print(scansion_result)

perform_scansion()