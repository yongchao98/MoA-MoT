def calculate_eoc_health():
    """
    Calculates and explains the theoretical maximum health of the Eye of Cthulhu
    in Terraria's Expert Mode with an infinite number of players.
    """

    # According to the official Terraria Wiki, which references the game's source code,
    # the health of a boss in Expert Mode is determined by a scaling formula.
    
    # 1. Define the base parameters for the calculation.
    base_health_eoc = 3640  # Base health of Eye of Cthulhu in Expert Mode.
    player_health_scaler = 0.35  # Health increases by 35% for each player after the first.
    num_players_symbol = "N" # We use 'N' to represent the number of players.

    # 2. Explain the formula and its components.
    print("The formula for boss health in Expert Mode multiplayer is:")
    print("Total Health = Base Health * (1 + (Number of Players - 1) * Player Health Scaler)")
    print("\nFor the Eye of Cthulhu, the specific values are:")
    print(f"Base Health = {base_health_eoc}")
    print(f"Player Health Scaler = {player_health_scaler}")
    print(f"Number of Players = {num_players_symbol}")

    # 3. Present the final equation with the numbers substituted.
    print("\nPlugging these numbers into the formula, we get the final equation:")
    print(f"Total Health = {base_health_eoc} * (1 + ({num_players_symbol} - 1) * {player_health_scaler})")

    # 4. Analyze the equation as the number of players approaches infinity.
    print("\nIn this case, the number of players (N) is considered to be infinite.")
    print("As 'N' approaches infinity, the term '(N - 1)' also approaches infinity.")
    print("When this infinitely large number is multiplied by the positive constants (0.35 and 3640), the result is also infinity.")
    print("\nTherefore, the theoretical maximum health of the Eye of Cthulhu with an infinite number of players is infinite.")

if __name__ == "__main__":
    calculate_eoc_health()