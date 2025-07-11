import sys

def solve():
    """
    Reads a single word and applies the abbreviation rule if it's too long.
    A word is too long if its length is strictly greater than 8 characters.
    """
    try:
        # Read a single line from standard input, which represents the word.
        # .strip() removes leading/trailing whitespace, including the newline character.
        word = sys.stdin.readline().strip()

        # We do nothing for empty input
        if not word:
            return

        # Get the length of the word.
        length = len(word)

        # The rule applies if the length is strictly more than 8.
        if length > 8:
            # The abbreviation consists of the first letter...
            first_letter = word[0]
            # ...the last letter...
            last_letter = word[-1]
            # ...and the count of letters between them.
            # This is the "number in the final equation" as requested.
            count_in_between = length - 2
            
            # Construct the final abbreviated string and print it.
            # For example, "localization" -> "l" + "10" + "n"
            print(f"{first_letter}{count_in_between}{last_letter}")
        else:
            # If the word is not too long, print it as is.
            print(word)

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)

solve()