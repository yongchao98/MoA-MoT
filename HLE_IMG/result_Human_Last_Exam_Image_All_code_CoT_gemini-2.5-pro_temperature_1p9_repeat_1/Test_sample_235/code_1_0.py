def solve():
    """
    This function calculates the total syllables in the erasure poem to identify its metric pattern.
    """
    # The words of the poem are interpreted to include the title, page number, and collage elements.
    # The number "35" is spelled out as "thirty-five".
    poem_words = ["The", "first", "step", "thirty-five", "rules", "and", "lines", "an", "intricate", "spider's", "web", "work"]

    # A dictionary mapping each word to its syllable count.
    syllable_counts = {
        "The": 1,
        "first": 1,
        "step": 1,
        "thirty-five": 3,
        "rules": 1,
        "and": 1,
        "lines": 1,
        "an": 1,
        "intricate": 3,
        "spider's": 2,
        "web": 1,
        "work": 1
    }

    total_syllables = 0
    calculation_string_parts = []

    print("Calculating the total number of syllables in the poem:")
    
    for word in poem_words:
        count = syllable_counts.get(word, 0)
        print(f"- '{word}': {count} syllable(s)")
        total_syllables += count
        calculation_string_parts.append(str(count))
    
    # Joining the numbers with ' + ' for the final equation.
    calculation_string = " + ".join(calculation_string_parts)

    print("\nFinal calculation:")
    print(f"{calculation_string} = {total_syllables}")
    print("\nThe poem has a total of 17 syllables. A single sentence of 17 syllables is known as an 'American Sentence'.")

solve()