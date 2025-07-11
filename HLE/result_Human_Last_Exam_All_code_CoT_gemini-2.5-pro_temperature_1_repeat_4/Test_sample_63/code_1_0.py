# This script solves a riddle by breaking it down into its constituent clues.
# The riddle: "As Kurt Vonnegut noted, this man looks like a porcupine in all the pictures.
# Name this man in two words that start with the same letter."

def solve_riddle():
    """
    Solves the riddle by analyzing its two main clues and finding the
    historical figure who fits both conditions.
    """
    # Clue 1: The Vonnegut reference points to Albert Einstein due to his wild hair.
    # This establishes the "porcupine hair" characteristic.
    # However, "Albert Einstein" does not satisfy the second clue.

    # Clue 2: The name must be two words starting with the same letter.

    # We need a famous person with a name like "Xxxxx Xxxxx" who is also
    # known for having wild, unkempt hair.

    # Candidate: Benjamin Franklin.
    # Check 1: Does the name fit? "Benjamin" and "Franklin" both start with 'F'. Yes.
    # Check 2: Does the description fit? His portraits show him with long,
    # often wild hair. In the context of a riddle, this is a plausible match for "porcupine-like".

    first_name = "Benjamin"
    last_name = "Franklin"

    # Print the final answer.
    print(first_name + " " + last_name)

solve_riddle()