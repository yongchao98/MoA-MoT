import sys

def solve():
    """
    This function reads words and applies the abbreviation rule.
    """
    # The problem description implies reading multiple words, one per line.
    # The C code, however, processes only one line.
    # To be robust, this solution will handle multiple lines of input
    # until there is no more input.
    
    # Let's read the number of words first, as is common in competitive programming.
    try:
        n_str = sys.stdin.readline()
        if not n_str:
            return
        n = int(n_str.strip())
    except (ValueError, IndexError):
        # If the first line is not a number, treat it as the word itself.
        # This matches the behavior of the C code which processes one line.
        process_word(n_str)
        return

    for _ in range(n):
        word = sys.stdin.readline().strip()
        if word:
            process_word(word)

def process_word(word):
    """
    Processes a single word and prints either the word or its abbreviation.
    The rule is to abbreviate if the length is strictly greater than 8.
    """
    length = len(word)
    
    # The problem states to abbreviate if length is "strictly more than eight characters".
    if length > 8:
        first_letter = word[0]
        last_letter = word[-1]
        
        # The number of letters between the first and the last.
        middle_count = length - 2
        
        # The final output is the first letter, the count, and the last letter.
        # The prompt asks to "output each number in the final equation".
        # We will print the components that form the final abbreviated word.
        # The number `middle_count` is printed as a single decimal value.
        print(f"{first_letter}{middle_count}{last_letter}")
    else:
        # If the word is not too long, print it as is.
        print(word)

# In a typical competitive programming setup, the main logic is called directly.
# The C code provided processes only one line, so we will write a version
# that also processes just one line to be analogous.

def solve_single_word():
    """
    This function solves the problem for a single word input,
    matching the behavior of the provided C code.
    """
    word = sys.stdin.readline().strip()
    if not word:
        return
        
    length = len(word)
    
    if length > 8:
        first_letter = word[0]
        last_letter = word[-1]
        middle_count = length - 2
        
        # Printing the final abbreviation string.
        # The number is middle_count.
        print(f"{first_letter}{middle_count}{last_letter}")
    else:
        print(word)

# Execute the single-word solution.
solve_single_word()