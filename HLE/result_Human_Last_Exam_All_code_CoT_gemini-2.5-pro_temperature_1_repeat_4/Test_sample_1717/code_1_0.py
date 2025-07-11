def solve_terraria_health_limit():
    """
    Calculates and explains the theoretical maximum health of the Eye of Cthulhu
    in Terraria's Expert Mode with an infinite number of players.
    """
    # Step 1: Define the base health for the Eye of Cthulhu in Expert Mode.
    base_health_expert = 3640

    # Step 2: Define the health increase factor for each additional player.
    # In Expert Mode, boss health increases by 35% for each player after the first.
    health_increase_factor = 0.35

    # The formula for health (H) based on the number of players (N) is:
    # H = BaseHealth * (1 + HealthIncreaseFactor * (N - 1))
    # As N approaches infinity, H also approaches infinity.

    # Step 3: Identify the game's technical limitation.
    # Terraria's source code stores NPC health in a standard 32-bit signed integer (int).
    # This data type has a fixed maximum value.
    theoretical_max_health = 2147483647

    # Step 4: Print the explanation and the final answer.
    print("In Terraria's Expert Mode, the health of the Eye of Cthulhu scales with the number of players.")
    print("\nThe formula is: Health = Base Health * (1 + Per-Player Bonus * (Number of Players - 1))")
    print("\nUsing the specific numbers for the Eye of Cthulhu:")
    print(f"Health = {base_health_expert} * (1 + {health_increase_factor} * (Number of Players - 1))")
    
    print("\nWith an infinite number of players, this formula would yield infinite health.")
    print("However, the game's code stores health as a 32-bit signed integer, which has a maximum value.")
    print("This value acts as the true theoretical health limit in the game engine.")
    
    print("\nTherefore, the theoretical maximum health is:")
    print(f"{theoretical_max_health}")

solve_terraria_health_limit()