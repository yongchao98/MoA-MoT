import sys

def solve():
    """
    This function solves the too long words problem.
    It reads words from standard input and applies the abbreviation rule.
    """
    try:
        # According to standard competitive programming format for this problem,
        # the first line is an integer n, the number of words.
        # This part is omitted in the C code but is typical for the problem.
        # If we are to process just one word, this loop would not be necessary.
        # The C code is written to process one word only. Here we follow the C
        # code's behavior and process a single line of input.
        
        word = sys.stdin.readline().strip()
        
        if not word:
            return

        n = len(word)
        if n > 8:
            # Abbreviate the word
            first_char = word[0]
            last_char = word[-1]
            middle_count = n - 2
            # The following print statement outputs the required format,
            # including the number in the middle of the new string.
            print(f"{first_char}{middle_count}{last_char}")
        else:
            # The word is not too long, print it as is
            print(word)

    except (IOError, ValueError) as e:
        # Handle potential input errors gracefully
        print(f"An error occurred: {e}", file=sys.stderr)

# We can call the function to execute the logic.
# For a competitive programming environment, you would typically just have the code run directly.
# Let's assume we are running this as a script for a single word.
solve()
