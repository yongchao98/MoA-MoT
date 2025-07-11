def solve_terraria_health_cap():
    """
    Calculates and explains the theoretical maximum health of the Eye of Cthulhu
    in Expert Mode with an infinite number of players, based on Terraria's source code.
    """

    # Step 1: Define the base health for single-player Expert Mode.
    expert_base_health = 3640

    # Step 2: Define the health scaling factor per additional player.
    per_player_scaling_factor = 0.35

    # Step 3: Explain the health calculation formula.
    print("The health of the Eye of Cthulhu in Expert Mode is calculated based on the number of players (N).")
    print(f"Base Health (1 Player): {expert_base_health}")
    print(f"Health increase per additional player: {per_player_scaling_factor * 100}%")
    print(f"Full Formula: Total Health = {expert_base_health} * (1 + {per_player_scaling_factor} * (N - 1))\n")

    # Step 4: Explain the implication of infinite players and the technical limitation.
    print("As the number of players (N) approaches infinity, the result of this formula also approaches infinity.")
    print("However, Terraria's source code stores NPC health as a 32-bit signed integer.")
    print("This data type has a hard-coded maximum value, which acts as the game's theoretical health cap.\n")

    # Step 5: State the final answer, which is the maximum value for a 32-bit signed integer.
    # This value is 2^31 - 1.
    max_health_cap = 2147483647

    # Final Answer: Present the equation representing the limit.
    print("Final Equation:")
    print(f"Theoretical Maximum Health = Limit (32-bit Integer) = {max_health_cap}")


solve_terraria_health_cap()
<<<2147483647>>>