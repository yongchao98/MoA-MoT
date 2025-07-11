def solve_task():
    """
    This function generates and prints a phrase according to the user's specifications.
    """
    # 1. Define the words of the phrase and their corresponding syllable counts.
    # The phrase is "Unchecked Cell Division In Luys Body Is Ferocious".
    words = ["Unchecked", "Cell", "Division", "In", "Luys", "Body", "Is", "Ferocious"]
    syllables = [2, 1, 3, 1, 1, 2, 1, 3]

    # 2. Combine the words to form the final, capitalized phrase.
    # The words are already capitalized and contain no hyphens or apostrophes.
    final_phrase = " ".join(words)

    # 3. Create a string representing the "equation" of the syllable counts.
    equation_parts = [str(s) for s in syllables]
    total_syllables = sum(syllables)
    final_equation = f"{' + '.join(equation_parts)} = {total_syllables}"

    # 4. Print the final phrase and the equation.
    print(final_phrase)
    print("Syllable count equation:")
    print(final_equation)

solve_task()