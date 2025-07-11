def calculate_terraria_health_limit():
    """
    Calculates the theoretical maximum health of the Eye of Cthulhu in Expert Mode
    with an infinite number of players, based on data type limitations in the game's code.
    """

    # Step 1: Define the base parameters for the Eye of Cthulhu's health calculation.
    base_health_classic = 2800
    expert_mode_multiplier = 1.3
    health_per_player_scaling_factor = 0.35
    base_health_expert = int(base_health_classic * expert_mode_multiplier)

    # Step 2: Explain the scaling formula and the concept of the theoretical limit.
    print("In Terraria's Expert Mode, a boss's health scales with the number of players.")
    print(f"The Eye of Cthulhu's base health in Expert Mode for one player is {base_health_expert}.")
    print(f"For each additional player, the total health is increased by a factor of {health_per_player_scaling_factor} (35%).")
    print(f"The formula is: Total Health = {base_health_expert} * (1 + {health_per_player_scaling_factor} * (PlayerCount - 1))\n")
    print("With an infinite number of players, this formula would mathematically approach infinity.")
    print("However, the game's source code stores NPC health in a 32-bit signed integer.")
    print("This data type has a hard-coded maximum value, which represents the true theoretical maximum health for any entity in the game.\n")

    # Step 3: Calculate and display the final equation for this maximum value.
    # The maximum value for a 32-bit signed integer is 2 to the power of 31, minus 1.
    power = 31
    base = 2
    subtract = 1
    theoretical_max_health = (base ** power) - subtract

    print("The final equation for the theoretical maximum health is determined by this data type limit:")
    print(f"Theoretical Max Health = {base} ** {power} - {subtract} = {theoretical_max_health:,}")


calculate_terraria_health_limit()