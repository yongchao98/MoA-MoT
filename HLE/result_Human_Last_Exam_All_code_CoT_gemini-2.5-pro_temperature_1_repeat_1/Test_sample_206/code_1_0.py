def solve_mtg_puzzle():
    """
    This function determines the optimal attack and resulting casualties for the described MTG scenario.
    """

    # Based on the analysis, the optimal play is to attack with only one creature.
    attacking_creatures = [4]

    # This attack forces a trade with the opponent's only flying creature.
    # Both creatures will die in combat.
    creatures_that_die = [4, 8]

    # The problem requires the numbers to be in increasing order.
    attacking_creatures.sort()
    creatures_that_die.sort()

    # Format the output as specified: (attacker1), (attacker2); (dying1), (dying2)
    # The join function handles the cases of single or multiple creatures in each list.
    attack_string = ", ".join(f"({c})" for c in attacking_creatures)
    die_string = ", ".join(f"({c})" for c in creatures_that_die)

    # Print the final formatted answer. The print statement outputs each required number.
    print(f"{attack_string}; {die_string}")

solve_mtg_puzzle()