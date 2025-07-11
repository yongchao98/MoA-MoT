import math

def solve_terraria_health():
    """
    Calculates and explains the theoretical maximum health of the Eye of Cthulhu
    in Terraria's Expert Mode with an infinite number of players.
    """

    # Step 1 & 2: Define constants based on Terraria's mechanics.
    # In Expert Mode, the Eye of Cthulhu's base health for one player is 3640.
    base_health_expert = 3640

    # In Expert Mode, health increases by 35% of the base for each player after the first.
    per_player_increase_multiplier = 0.35

    # Step 3: Present the formula.
    print("The health scaling formula for a boss in Expert Mode is:")
    print("Total Health = BaseHealth * (1 + (NumberOfPlayers - 1) * PlayerMultiplier)")
    print("\nFor the Eye of Cthulhu, we plug in the specific numbers:")

    # As requested, printing the final equation with each number.
    print(f"Total Health = {base_health_expert} * (1 + (NumberOfPlayers - 1) * {per_player_increase_multiplier})")

    # Step 4: Analyze for an infinite number of players.
    print("\nThis formula shows that the Total Health increases linearly with the 'NumberOfPlayers'.")
    print("As the 'NumberOfPlayers' approaches infinity, the 'Total Health' also increases without any upper limit.")
    print("\nTherefore, the theoretical maximum health is infinite.")

solve_terraria_health()