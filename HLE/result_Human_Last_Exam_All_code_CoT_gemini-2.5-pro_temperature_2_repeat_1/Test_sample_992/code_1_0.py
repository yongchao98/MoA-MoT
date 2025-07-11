def solve_character_riddle():
    """
    This function analyzes the structure of the Chinese character 森 (sēn)
    and prints an equation representing its stroke count.
    """
    # The radical 木 (mù, meaning tree) is the building block.
    # We define its properties.
    radical_strokes = 4
    num_radicals_in_character = 3

    # The character 森 (sēn, meaning forest) is made of three 木 radicals.
    total_strokes = radical_strokes * num_radicals_in_character

    # Build the equation string "4 + 4 + 4"
    equation_parts = [str(radical_strokes)] * num_radicals_in_character
    equation_sum_part = " + ".join(equation_parts)

    # Print the explanation and the final equation.
    print(f"The riddle refers to the character 森 (sēn), meaning forest.")
    print(f"It is composed of {num_radicals_in_character} 木 (mù) radicals, each with {radical_strokes} strokes.")
    print("The equation for the total number of strokes is:")
    print(f"{equation_sum_part} = {total_strokes}")

solve_character_riddle()