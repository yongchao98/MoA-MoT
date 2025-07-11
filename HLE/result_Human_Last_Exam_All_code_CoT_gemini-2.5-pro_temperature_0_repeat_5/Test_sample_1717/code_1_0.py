def calculate_max_health():
    """
    Calculates and explains the theoretical maximum health of the Eye of Cthulhu
    in Terraria Expert Mode with an infinite number of players.
    """

    # Step 1: Define the base parameters for the calculation.
    # In Expert Mode, the Eye of Cthulhu's base health is 3640.
    # Health increases by 35% for each player after the first.
    base_health = 3640
    per_player_increase = 0.35

    # Step 2: Explain the health scaling formula.
    print("The health scaling formula for a boss in Expert Mode multiplayer is:")
    print(f"Total Health = Base Health * (1 + Per-Player Increase * (Number of Players - 1))")
    print("Substituting the values for the Eye of Cthulhu:")
    print(f"Total Health = {base_health} * (1 + {per_player_increase} * (Number of Players - 1))\n")

    # Step 3: Explain the implication of infinite players and the technical limitation.
    print("As you can see, health scales linearly with the number of players.")
    print("Mathematically, with infinite players, the health would also be infinite.")
    print("However, the game's source code uses a 32-bit signed integer to store an NPC's health.")
    print("This data type has a hard-coded maximum value, which acts as the true ceiling for health in the game engine.\n")

    # Step 4: State the final maximum health based on the 32-bit integer limit.
    # The maximum value for a 32-bit signed integer is 2^31 - 1.
    max_int_32 = 2147483647

    print("The final equation for the theoretical maximum health is therefore the limit of this data type.")
    print(f"Maximum Health = {max_int_32}")


if __name__ == "__main__":
    calculate_max_health()