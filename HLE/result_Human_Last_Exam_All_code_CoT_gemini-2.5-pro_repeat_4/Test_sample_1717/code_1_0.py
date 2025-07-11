def calculate_eoc_health():
    """
    Calculates and explains the theoretical maximum health of the Eye of Cthulhu
    in Terraria's Expert Mode with an infinite number of players.
    """

    # Step 1: Define the base parameters for the Eye of Cthulhu in Expert Mode.
    base_health = 3640
    health_increase_per_player = 0.35  # This represents a 35% increase.

    # Step 2: Explain the health scaling formula.
    # The formula is: Total Health = BaseHealth * (1 + HealthIncrease * (NumberOfPlayers - 1))
    # Let's represent the number of players with the variable 'N'.
    print("The health scaling formula for an Expert Mode boss is:")
    print("Total Health = Base Health * (1 + Health Increase Per Player * (N - 1))\n")

    print("For the Eye of Cthulhu, the specific formula is:")
    print(f"Total Health = {base_health} * (1 + {health_increase_per_player} * (N - 1))\n")

    # Step 3: Simplify the formula into a linear equation (y = mx + b form).
    # Total Health = 3640 * (1 + 0.35*N - 0.35)
    # Total Health = 3640 * (0.65 + 0.35*N)
    # Total Health = (3640 * 0.65) + (3640 * 0.35 * N)
    constant_part = base_health * (1 - health_increase_per_player)
    scaling_part = base_health * health_increase_per_player
    
    print("This equation can be simplified to a linear equation of the form: y = mx + b")
    print("Where 'y' is the Total Health and 'x' is the number of players 'N'.\n")

    print("The simplified final equation is:")
    # The numbers in the final equation are explicitly printed as requested.
    print(f"Total Health = {int(constant_part)} + {int(scaling_part)} * N\n")

    # Step 4: Explain the conclusion for an infinite number of players.
    print("Conclusion:")
    print("As you can see from the final equation, the boss's health increases linearly with the number of players (N).")
    print("There is no upper limit or cap to this scaling in the game's code.")
    print("Therefore, as the number of players 'N' approaches infinity, the Total Health also approaches infinity.")

calculate_eoc_health()