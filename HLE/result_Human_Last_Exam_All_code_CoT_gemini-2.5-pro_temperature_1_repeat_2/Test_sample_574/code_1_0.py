def solve_puzzle():
    """
    This script solves the word puzzle by identifying the components
    of a Culture series ship name based on the given clues.
    """

    # Clue 1: Sleeveless garments that drape over the back and shoulders -> ROBES
    # Word formed from the letters of ROBES -> SOBER
    part1 = "SOBER"

    # Clue 2: Experienced and trusted individuals who guide and advise others -> COUNSEL
    # Word formed from the letters of COUNSEL -> COUNSEL
    part2 = "COUNSEL"

    # Combine the parts to form the final ship name
    ship_name = f"{part1} {part2}"

    # Print the equation showing how the final name is constructed
    print(f"The first word, derived from a sleeveless garment (ROBES), is: {part1}")
    print(f"The second word, from a term for trusted advisors (COUNSEL), is: {part2}")
    print("\nCombining them gives the Culture ship name:")
    print(f'"{part1}" + "{part2}" = "{ship_name}"')

solve_puzzle()