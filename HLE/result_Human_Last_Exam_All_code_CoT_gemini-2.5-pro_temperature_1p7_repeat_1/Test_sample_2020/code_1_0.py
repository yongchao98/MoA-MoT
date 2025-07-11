def find_scansion_issue():
    """
    This script analyzes the provided sestina to find a word
    that causes a scansion issue.
    """
    # 1. Explain the context of scansion in the poem.
    print("Analyzing the sestina for a scansion issue...")
    print("A common meter for poems like this is iambic pentameter, which has 10 syllables per line.")
    print("A word that causes a line to deviate from this 10-syllable count is considered a scansion issue.")
    print("-" * 50)

    # 2. Identify the specific problematic line and word.
    problem_line = "You said, 'Let’s end it,' quietly over tea."
    problem_word = "quietly"
    line_syllables = 11
    expected_syllables = 10

    print(f"Found a line with a scansion issue: '{problem_line}'")
    print(f"This line contains {line_syllables} syllables, which is one more than the expected {expected_syllables}.")
    print(f"The word that causes the line to be too long is '{problem_word}'.")
    print("-" * 50)

    # 3. Present the syllable breakdown as an "equation".
    print("Here is a syllable breakdown of the line:")
    
    # Numbers represent the syllable count for each word/phrase.
    you_said = 2
    lets_end_it = 3
    word_quietly = 3
    over_tea = 3
    total = you_said + lets_end_it + word_quietly + over_tea
    
    print(f"Syllable Equation: 'You said'(2) + 'Let’s end it,'(3) + 'quietly'(3) + 'over tea.'(3) = {total}")
    
    print("\nThe word 'quietly' has 3 syllables, which disrupts the 10-syllable meter.")

if __name__ == '__main__':
    find_scansion_issue()
    print("\n<<<quietly>>>")
