import collections

def find_poetic_form():
    """
    Identifies the poetic form by analyzing the syllable count of the given lines.
    """
    print("Analyzing the structure of the poem's lines...\n")

    # According to the prompt, the third line is from the image, 
    # and the final line is also given.
    # Line 3: "ghostly velum forms like a dance"
    # Line 4: "nacreous wavers"

    # Step 1 & 2: Create a breakdown of each line with its syllable count per word.
    # Using an OrderedDict to maintain the order of words for printing.
    line_3_word_syllables = collections.OrderedDict([
        ("ghost-ly", 2),
        ("vel-um", 2),
        ("forms", 1),
        ("like", 1),
        ("a", 1),
        ("dance", 1)
    ])
    
    line_4_word_syllables = collections.OrderedDict([
        ("na-cre-ous", 3),
        ("wa-vers", 2)
    ])

    # Step 3: Calculate and display the total syllables for Line 3.
    syllables_line_3 = sum(line_3_word_syllables.values())
    calculation_string_3 = " + ".join(map(str, line_3_word_syllables.values()))
    print("Syllable count for Line 3 ('ghostly velum forms like a dance'):")
    # As requested, printing each number in the final equation
    print(f"The sum of syllables is: {calculation_string_3} = {syllables_line_3}")
    print("-" * 40)

    # Step 3: Calculate and display the total syllables for Line 4.
    syllables_line_4 = sum(line_4_word_syllables.values())
    calculation_string_4 = " + ".join(map(str, line_4_word_syllables.values()))
    print("Syllable count for Line 4 ('nacreous wavers'):")
    # As requested, printing each number in the final equation
    print(f"The sum of syllables is: {calculation_string_4} = {syllables_line_4}")
    print("-" * 40)

    # Step 4 & 5: Analyze the pattern and conclude.
    fibonacci_sequence = [0, 1, 1, 2, 3, 5, 8, 13, 21]
    print(f"The syllable counts for the third and fourth lines are {syllables_line_3} and {syllables_line_4}.")
    print(f"Let's compare this to the Fibonacci sequence: {fibonacci_sequence}...")
    
    if syllables_line_3 in fibonacci_sequence and syllables_line_4 in fibonacci_sequence:
        print(f"\nBoth syllable counts, {syllables_line_4} and {syllables_line_3}, are consecutive numbers in the Fibonacci sequence.")
        print("A poem where the syllable count of the lines follows this sequence is called a 'Fibonacci poem' or 'Fib'.")
        print("\nConclusion: The text adheres to the form of a Fibonacci poem.")
    else:
        print("\nThe syllable counts do not correspond to a standard poetic form based on this analysis.")

find_poetic_form()