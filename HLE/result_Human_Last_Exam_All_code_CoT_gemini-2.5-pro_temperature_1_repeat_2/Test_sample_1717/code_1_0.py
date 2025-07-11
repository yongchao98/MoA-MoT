import math

def solve_terraria_health_limit():
    """
    Calculates and explains the theoretical maximum health for the Eye of Cthulhu
    in Terraria's Expert Mode with an infinite number of players.
    """
    # Step 1: Define the constants from the Terraria health scaling formula.
    base_health = 3640  # Base health for the Eye of Cthulhu in Expert Mode.
    per_player_increase = 0.35  # Health increases by 35% for each additional player.

    # Step 2: Define the technical limitation from the game's source code.
    # In Terraria (C#), NPC health is stored as a 32-bit signed integer (System.Int32).
    # This data type has a maximum possible value.
    max_health_limit = 2**31 - 1  # This is 2,147,483,647.

    # Step 3: Explain the logic.
    # The formula for health is: Total Health = base_health * (1 + per_player_increase * (NumberOfPlayers - 1))
    # As the number of players approaches infinity, the result of this formula also approaches infinity.
    # However, the health value is capped by the maximum value of the integer data type used to store it.

    print("Terraria Boss Health Calculation Analysis")
    print("-" * 40)
    print("The health scaling formula for a boss in Expert Mode is:")
    print("Final Health = Base Health * (1 + Player Increase * (Number of Players - 1))")
    print("\nBreaking down the formula for the Eye of Cthulhu:")
    # The final prompt requires printing each number in the final equation.
    print(f"Final Health = {base_health} * (1 + {per_player_increase} * (Number of Players - 1))")

    print("\nWith an infinite number of players, the result of this formula would mathematically be infinite.")
    print("However, the game's code stores health in a 32-bit signed integer, which has a hard technical limit.")
    print("This limit is the highest number the variable can hold before it overflows.")

    print("\nTheoretical Maximum Health:")
    # Step 4: Print the final answer.
    print(f"The theoretical maximum health is capped by the 32-bit integer limit, which is {max_health_limit}.")

solve_terraria_health_limit()