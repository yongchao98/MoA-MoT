import re

def solve_reading_time_puzzle():
    """
    Solves a psycholinguistic puzzle to find the word position
    of expected reading time elevation.
    """
    passage = "Experience will be a decisive factor in this battle. The ships will outmaneuver the enemy fleet because the old man the boats."
    critical_word = "man"

    # Normalize and split the passage into words.
    # Remove punctuation and convert to lowercase.
    cleaned_passage = re.sub(r'[^\w\s]', '', passage).lower()
    words = cleaned_passage.split()

    # Find the 1-based index (ordinal position) of the critical word.
    try:
        # The list is 0-indexed, so we add 1 for the ordinal position.
        ordinal_position = words.index(critical_word) + 1
    except ValueError:
        print(f"The critical word '{critical_word}' was not found in the passage.")
        return

    # A simple dictionary to convert numbers to words as requested by the prompt.
    num_to_word = {
        1: "one", 2: "two", 3: "three", 4: "four", 5: "five",
        6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten",
        11: "eleven", 12: "twelve", 13: "thirteen", 14: "fourteen",
        15: "fifteen", 16: "sixteen", 17: "seventeen", 18: "eighteen",
        19: "nineteen", 20: "twenty", 21: "twenty-one", 22: "twenty-two",
        23: "twenty-three", 24: "twenty-four", 25: "twenty-five"
        # Extend this dict if needed for other problems
    }

    # "output each number in the final equation" - Here we output the calculated number.
    print(f"The word causing the slowdown is '{critical_word}'.")
    print(f"The calculated ordinal word position is: {ordinal_position}")

    # Convert the number to a word for the final answer.
    final_answer = num_to_word.get(ordinal_position, str(ordinal_position))

    # The prompt requests the final answer formatted in a specific way.
    print(f"The answer as a single word in lowercase letters is: {final_answer}")

solve_reading_time_puzzle()