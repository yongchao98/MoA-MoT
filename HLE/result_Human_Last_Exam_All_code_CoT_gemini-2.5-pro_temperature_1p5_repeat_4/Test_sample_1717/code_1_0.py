def solve_terraria_health_limit():
    """
    Calculates and explains the theoretical maximum health for the Eye of Cthulhu
    in Terraria's Expert Mode, considering game engine limitations.
    """

    # Parameters from Terraria's source for Eye of Cthulhu (Expert Mode)
    base_health = 3640
    player_scaling_factor = 0.35

    # The health increases by a fixed amount for each player after the first.
    health_increase_per_player = int(base_health * player_scaling_factor)

    # --- Step 1: Explain and print the health scaling equation ---
    print("In Terraria's Expert Mode, boss health scales with the number of players, 'n'.")
    print("The formula is: Total Health = Base Health + (n - 1) * Health Increase Per Player")
    print("\nFor the Eye of Cthulhu, the specific equation is:")
    # The user requested to output each number in the final equation.
    # The equation for 'n' players is the most relevant one to show.
    print(f"Total Health = {base_health} + (n - 1) * {health_increase_per_player}")
    print("-" * 30)

    # --- Step 2: Explain the implication of "infinite players" and the code limit ---
    print("With the above formula, an infinite number of players would lead to infinite health.")
    print("\nHowever, the question refers to the 'source code'. In Terraria's code, an NPC's health is stored as a 32-bit signed integer.")
    print("This data type has a hard-coded maximum value, which imposes the actual theoretical limit within the game engine.")
    print("-" * 30)

    # --- Step 3: State the final answer ---
    # The maximum value for a 32-bit signed integer (System.Int32 in C#)
    max_health_value = 2147483647

    print("The final equation for the theoretical maximum health is therefore determined by this limit:")
    print(f"Theoretical Maximum Health = {max_health_value}")


solve_terraria_health_limit()