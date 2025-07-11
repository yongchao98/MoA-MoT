def solve_task():
    """
    This function generates a fourteen-syllable phrase according to the user's request
    and provides a syllable breakdown to meet all constraints.
    """
    
    # Define the words of the phrase and their syllable counts.
    words = ["Apoptosis", "In", "The", "Luys", "Results", "In", "Growth", "Ferocious"]
    syllables = [4, 1, 1, 1, 2, 1, 1, 3]

    # Construct the final phrase with each word capitalized.
    final_phrase = " ".join(words)
    
    # Construct the syllable "equation" string to demonstrate the count.
    # This fulfills the "output each number in the final equation" requirement.
    total_syllables = sum(syllables)
    equation_string_parts = [str(s) for s in syllables]
    equation_string = " + ".join(equation_string_parts) + f" = {total_syllables}"

    # Print the final phrase and the syllable equation.
    print(final_phrase)
    print("\nSyllable Equation:")
    print(equation_string)

solve_task()