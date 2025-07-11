def explain_eoc_max_health():
    """
    Explains and calculates the theoretical maximum health of the Eye of Cthulhu
    in Terraria's Expert Mode with an infinite number of players.
    """

    # According to the Terraria source code and wiki, here are the base values for Expert Mode.
    base_health = 4275
    health_multiplier_per_extra_player = 0.35  # This represents a 35% increase.

    print("To find the theoretical maximum health of the Eye of Cthulhu in Expert Mode, we need to use its health scaling formula.")
    print("-" * 50)
    print("The formula for boss health in Expert Mode multiplayer is:")
    print("Total Health = Base Health * (1 + Scaling Factor * (Number of Players - 1))\n")

    print("For the Eye of Cthulhu, the specific numbers are:")
    print(f"Base Health: {base_health}")
    print(f"Scaling Factor: {health_multiplier_per_extra_player} (or 35%) per player after the first.")
    print("-" * 50)

    print("The final equation is:")
    print(f"Total Health = {base_health} * (1 + {health_multiplier_per_extra_player} * (Number of Players - 1))\n")

    print("The question asks for the health with an infinite number of players.")
    print("As the 'Number of Players' in the equation approaches infinity, the 'Total Health' also increases without any limit.")
    print("Therefore, the theoretical maximum health is not a finite number.")
    print("\nConclusion: The theoretical maximum health is infinite.")

explain_eoc_max_health()