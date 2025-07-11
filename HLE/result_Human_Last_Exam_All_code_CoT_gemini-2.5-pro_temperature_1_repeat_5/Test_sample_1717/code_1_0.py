def calculate_eoc_health():
    """
    Calculates and explains the theoretical maximum health of the Eye of Cthulhu
    in Terraria's Expert Mode with an infinite number of players.
    """

    # Step 1: Define the base health for the Eye of Cthulhu in Expert Mode.
    # This value is from the game's source code for a single-player world.
    base_expert_health = 3640

    # Step 2: Define the multiplayer health scaling factor for Expert Mode.
    # For each player beyond the first, the boss's health multiplier increases by 0.35.
    per_player_increase = 0.35

    # Step 3: Construct the full health formula based on the number of players (N).
    # The formula is Health = BaseHealth * (1 + PerPlayerIncrease * (NumberOfPlayers - 1))
    print("The formula for the Eye of Cthulhu's health (H) in Expert Mode with N players is:")
    print(f"H(N) = {base_expert_health} * (1 + {per_player_increase} * (N - 1))")

    # To better analyze this, we can simplify the formula into a linear equation: y = mx + b
    # H(N) = 3640 * (1 + 0.35*N - 0.35)
    # H(N) = 3640 * (0.65 + 0.35*N)
    # H(N) = (3640 * 0.65) + (3640 * 0.35) * N
    # H(N) = 2366 + 1274 * N
    constant_term = base_expert_health * (1 - per_player_increase)
    linear_coefficient = base_expert_health * per_player_increase

    print("\nThis simplifies to the following linear equation:")
    print(f"H(N) = {int(constant_term)} + {int(linear_coefficient)} * N")

    # Step 4: Analyze the formula as N (number of players) approaches infinity.
    print("\nTo find the theoretical maximum health, we must find the limit of this function as N approaches infinity:")
    print(f"lim (N -> ∞) [ {int(constant_term)} + {int(linear_coefficient)} * N ]")
    
    print("\nAs the number of players (N) grows infinitely large, the term '1274 * N' also grows infinitely large.")
    print("Therefore, the total health increases without bound.")
    print("\nThe theoretical maximum health of the Eye of Cthulhu is infinity.")

calculate_eoc_health()