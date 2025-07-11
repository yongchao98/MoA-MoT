def solve_phrase_task():
    """
    This function generates and prints a phrase according to the user's specifications.
    It also shows the syllable calculation as an equation.
    """
    # Define the words of the phrase and their corresponding syllable counts.
    words = ["The", "Corpus", "Luysii's", "Tumor", "Growth", "Is", "So", "Atrocious"]
    syllables = [1, 2, 3, 2, 1, 1, 1, 3]

    # 1. Build and print the syllable equation.
    # This fulfills the requirement: "output each number in the final equation".
    equation_parts = [str(s) for s in syllables]
    total_syllables = sum(syllables)
    equation_string = " + ".join(equation_parts) + " = " + str(total_syllables)

    print("Syllable count breakdown for the 14-syllable phrase:")
    print(equation_string)
    print("-" * 30)

    # 2. Construct the final phrase with each word capitalized.
    final_phrase = " ".join(words)

    # 3. Print the final phrase.
    print("The final fourteen-syllable phrase is:")
    print(final_phrase)

# Execute the function to display the answer.
solve_phrase_task()
