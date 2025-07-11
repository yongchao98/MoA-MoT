def calculate_eoc_health():
    """
    Calculates the theoretical maximum health of the Eye of Cthulhu in Expert Mode.
    """
    # Step 1: Define the base health for the Eye of Cthulhu in single-player Expert Mode.
    # According to the official Terraria wiki, this value is 3640.
    base_expert_health = 3640

    # Step 2: Define the multiplayer health scaling factor.
    # In Expert Mode, boss health increases by 35% for each player after the first.
    per_player_increase = 0.35

    # Step 3: Present the formula for health based on the number of players, 'n'.
    # The formula is: Health = BaseHealth * (1 + PerPlayerIncrease * (n - 1))
    print("The formula for the Eye of Cthulhu's health in Expert Mode for 'n' players is:")
    print(f"Health(n) = {base_expert_health} * (1 + {per_player_increase} * (n - 1))")
    print("\nThis can be simplified to the linear equation:")
    
    # To simplify: Health(n) = 3640 * (1 + 0.35n - 0.35) = 3640 * (0.65 + 0.35n) = 2366 + 1274n
    constant_part = base_expert_health * (1 - per_player_increase)
    scaling_part = base_expert_health * per_player_increase
    print(f"Health(n) = {int(constant_part)} + {int(scaling_part)} * n")

    # Step 4: Analyze the formula for an infinite number of players.
    # As 'n' (the number of players) approaches infinity, the term '1274 * n' also approaches infinity.
    # A linear function with a positive slope increases without any upper bound.
    print("\nBecause the health scales linearly with the number of players, as the number of players approaches infinity, the health also increases without limit.")
    print("\nTherefore, the theoretical maximum health is infinite.")

calculate_eoc_health()