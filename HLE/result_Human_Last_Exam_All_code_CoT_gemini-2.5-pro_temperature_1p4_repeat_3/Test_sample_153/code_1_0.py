def solve_cold_war_puzzle():
    """
    Solves the puzzle by decoding the word "Кома" as a sum of atomic numbers.

    The steps are:
    1. The letters of "Кома" (К, О, М, А) are interpreted as the first letters
       of the Russian names for specific chemical elements.
    2. The chosen elements are Калий (Potassium), Кислород (Oxygen), 
       Марганец (Manganese), and Азот (Nitrogen).
    3. The atomic numbers of these elements are summed up.
    4. This sum corresponds to a Russian vehicle registration plate code, 
       which identifies the correct region.
    """

    # Dictionary mapping the Cyrillic letters to the chosen element and its atomic number.
    element_map = {
        'К': ('Potassium (Калий)', 19),
        'о': ('Oxygen (Кислород)', 8),
        'м': ('Manganese (Марганец)', 25),
        'а': ('Nitrogen (Азот)', 7)
    }

    # The clue word
    clue = "Кома"
    
    # Lists to store numbers and for the final equation string
    atomic_numbers = []
    equation_parts = []

    print("Decoding the clue 'Кома':")
    for char in clue:
        # We use lower() to match keys in the dictionary for consistency
        char_lower = char.lower()
        if char_lower in element_map:
            element_name, number = element_map[char_lower]
            print(f"  '{char}' -> {element_name}, Atomic Number: {number}")
            atomic_numbers.append(number)
            equation_parts.append(str(number))

    # Calculate the sum
    total = sum(atomic_numbers)
    
    # Build and print the final equation
    equation = " + ".join(equation_parts)
    print("\nCalculating the sum of the atomic numbers:")
    print(f"{equation} = {total}")

    print(f"\nThe resulting number, {total}, is the vehicle registration code for Perm Krai.")

solve_cold_war_puzzle()

<<<B>>>