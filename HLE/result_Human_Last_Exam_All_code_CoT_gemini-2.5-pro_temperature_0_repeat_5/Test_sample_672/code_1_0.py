import sys

# This script performs scansion on a line of poetry.
# The line is: "And in the letter, my cousin mentions a piece of advice"
# The scansion is determined by analyzing the natural stress pattern of the words.
# '/' represents a stressed syllable.
# 'x' represents an unstressed syllable.

def perform_scansion():
    """
    Analyzes the syllables and their stresses for the given line and prints the scansion pattern.
    """
    # Syllable-by-syllable stress analysis:
    # And: x (unstressed)
    # in: x (unstressed)
    # the: x (unstressed)
    # let-ter: /x (stressed, unstressed)
    # my: x (unstressed)
    # cou-sin: /x (stressed, unstressed)
    # men-tions: /x (stressed, unstressed)
    # a: x (unstressed)
    # piece: / (stressed)
    # of: x (unstressed)
    # ad-vice: x/ (unstressed, stressed)
    
    scansion_pattern = "xxx/xx/x/xx/xx/"
    
    print(scansion_pattern)

perform_scansion()