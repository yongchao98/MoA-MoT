def solve_riddle():
    """
    This function solves the historical riddle based on the provided clues.
    """
    # The clues point to Nikolai Rezanov.
    # 1. Voznesensky's rock opera "Juno and Avos" is about Rezanov's voyage.
    # 2. Rezanov's grave in Krasnoyarsk was lost and then memorialized in the late 1980s.
    last_name = "Rezanov"

    print(f"The man described by Andrei Voznesensky is Nikolai Rezanov.")
    print(f"His last name is: {last_name}")

    # As requested, here is a representation of the name as a 'final equation'
    # where each letter is assigned its value based on its position in the alphabet.
    print("\nA numerical equation for the name:")

    equation_parts = []
    # Create a string like "R(18) + E(5) + ..."
    for char in last_name.upper():
        # ord(char) gets the ASCII value. We subtract the value for 'A' and add 1.
        numerical_value = ord(char) - ord('A') + 1
        equation_parts.append(f"{char}({numerical_value})")

    # The final equation string showing each letter and its number
    final_equation = " + ".join(equation_parts)
    print(final_equation)

solve_riddle()