def solve_phrase_task():
    """
    This function solves the task by constructing a phrase based on a specific
    set of rules and verifying its syllable count with a printed equation.
    """
    
    # The words for the phrase are chosen to meet the prompt's constraints.
    words = ["A", "Malignant", "Growth", "And", "Awful", "Gloom", "Spell", "Luys'", "Body's", "Doom"]
    
    # The syllable count for each corresponding word.
    syllables = [1, 3, 1, 1, 2, 1, 1, 1, 2, 1]

    # Construct the final phrase by joining the words.
    # The prompt requires each word to be capitalized.
    final_phrase = " ".join(words)
    
    # Create the verification equation for the syllable count.
    # The prompt requires outputting each number in the final equation.
    equation_parts = [str(s) for s in syllables]
    total_syllables = sum(syllables)
    syllable_equation = " + ".join(equation_parts) + f" = {total_syllables}"

    # Print the final results as required.
    print("Generated Phrase:")
    print(final_phrase)
    print("\nSyllable Count Equation:")
    print(syllable_equation)

# Execute the function to display the solution.
solve_phrase_task()