import sys

def solve():
    """
    Reads a word, and if it's longer than 8 characters,
    abbreviates it and prints the result.
    Otherwise, prints the word as is.
    """
    # Read a single line from standard input and remove trailing whitespace
    word = sys.stdin.readline().strip()

    # Get the length of the word
    length = len(word)

    # Check if the word is "too long" (length > 8)
    if length > 8:
        # Get the first character
        first_char = word[0]
        # Get the last character
        last_char = word[-1]
        # Calculate the number of characters between the first and last
        middle_count = length - 2
        
        # Print the components of the abbreviation sequentially to form the final string.
        # This matches the behavior of the provided C code.
        print(first_char, end='')
        print(middle_count, end='')
        print(last_char)
    else:
        # If the word is not too long, print it as is
        print(word)

solve()