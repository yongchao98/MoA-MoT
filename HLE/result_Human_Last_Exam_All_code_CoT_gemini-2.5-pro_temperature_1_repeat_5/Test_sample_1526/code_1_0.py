import sys

def solve():
    """
    Reads words from standard input and applies the abbreviation rule.
    A word is abbreviated if its length is strictly greater than 8.
    The abbreviation consists of the first letter, the count of letters
    in between, and the last letter.
    """
    # The problem description is slightly ambiguous whether to process one word or many.
    # This implementation will process multiple words, one per line, until EOF.
    for line in sys.stdin:
        word = line.strip()
        if not word:
            continue

        n = len(word)

        # A word is "too long" if its length is strictly more than 8 characters.
        if n > 8:
            # The abbreviation is the first letter, the number of characters
            # between the first and last, and the last letter.
            first_char = word[0]
            last_char = word[-1]
            middle_count = n - 2
            
            # The problem asks to "output each number in the final equation!".
            # This is an unusual instruction. The C code and problem description
            # imply the output should be a single concatenated string like "l10n".
            # We will follow that concatenation logic.
            print(f"{first_char}{middle_count}{last_char}")
        else:
            # Words with length 8 or less are not abbreviated.
            print(word)

# To run this code, you would typically execute the script and provide
# input, for example:
# echo "localization" | python your_script_name.py
#
# Since the environment here does not execute the code, we will not call solve().
# If you wish to test it, uncomment the following line:
# solve()
# Example of what would be printed for "localization":
# l10n
# Example for "word":
# word
