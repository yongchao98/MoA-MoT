def solve_riddle():
    """
    Solves the riddle by identifying the word and the numerical clues within it.
    """
    # The word "X" that fits the historical and logical clues.
    # Being illiterate would be shameful for a Renaissance Pope.
    # The problem of illiteracy was widely "written" about in the 1960s (e.g., by UNESCO).
    answer_word = "ILLITERATE"

    # The riddle contains numerical hints ("Paul II") and the prompt asks for an equation.
    # The answer word itself contains Roman numerals.
    roman_numerals = {
        'I': 1,
        'L': 50,
        'C': 100,
        'D': 500,
        'M': 1000
    }

    # Find the numerals and their values within the word.
    found_numerals = []
    for char in answer_word:
        if char in roman_numerals:
            found_numerals.append((char, roman_numerals[char]))

    # Print the explanation and the "equation" as requested.
    print("The word 'ILLITERATE' contains the following Roman numerals:")

    equation_parts = []
    for char, value in found_numerals:
        equation_parts.append(f"{char} = {value}")

    # The riddle mentions Pope Paul II, and the word ILLITERATE has two I's and two L's.
    # We will display the values of each numeral found in the word.
    print("The 'equation' shows the value of each Roman numeral found in the word:")
    print("I = 1")
    print("L = 50")
    print("L = 50")
    print("I = 1")

    print("\nThe final answer is the word itself.")
    print(f"X = {answer_word}")


solve_riddle()
<<<ILLITERATE>>>