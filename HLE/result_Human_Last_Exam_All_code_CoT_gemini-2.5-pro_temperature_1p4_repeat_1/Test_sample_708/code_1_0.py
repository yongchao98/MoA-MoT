import re

def solve():
    """
    This script identifies the ordinal word position where an elevated reading time is expected
    in the given passage due to the successful resolution of a syntactic ambiguity.
    """
    passage = "Experience will be a decisive factor in this battle. The ships will outmaneuver the enemy fleet because the old man the boats."
    target_word = "man"

    # Step 1: Pre-process the passage to get a clean list of words.
    # We remove punctuation and split the string by spaces.
    words = re.sub(r'[^\w\s]', '', passage).split()

    # Step 2: Find the 1-based position of the target word.
    # Python's index() is 0-based, so we add 1.
    try:
        position = words.index(target_word) + 1
    except ValueError:
        print(f"The target word '{target_word}' was not found in the passage.")
        return

    # Step 3: Fulfill the instruction to show the numbers in the "final equation".
    # We interpret this as showing the calculation of the position.
    equation_str = " + ".join(["1"] * position)

    # Step 4: Convert the numerical position to an ordinal word.
    # We'll use a dictionary for common ordinals.
    ordinals = {
        1: 'first', 2: 'second', 3: 'third', 4: 'fourth', 5: 'fifth',
        6: 'sixth', 7: 'seventh', 8: 'eighth', 9: 'ninth', 10: 'tenth',
        11: 'eleventh', 12: 'twelfth', 13: 'thirteenth', 14: 'fourteenth',
        15: 'fifteenth', 16: 'sixteenth', 17: 'seventeenth', 18: 'eighteenth',
        19: 'nineteenth', 20: 'twentieth', 21: 'twenty-first', 22: 'twenty-second'
    }

    # A simple fallback for numbers not in our dictionary.
    if position in ordinals:
        ordinal_word = ordinals[position]
    elif str(position).endswith('1'):
        ordinal_word = f"{position}st"
    elif str(position).endswith('2'):
        ordinal_word = f"{position}nd"
    elif str(position).endswith('3'):
        ordinal_word = f"{position}rd"
    else:
        ordinal_word = f"{position}th"

    # Step 5: Print the results.
    print(f"Passage: \"{passage}\"")
    print(f"Tokenized words: {words}")
    print(f"The critical word where ambiguity is resolved: '{target_word}'")
    print(f"The ordinal position is calculated by counting the words.")
    print(f"Equation for position: {equation_str} = {position}")
    print(f"The final answer is the ordinal word for position {position}.")
    print(f"Final Answer: {ordinal_word}")

solve()