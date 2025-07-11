import math

def calculate_eoc_health():
    """
    Calculates and explains the theoretical maximum health of the Eye of Cthulhu
    in Terraria's Expert Mode with an infinite number of players.
    """
    # According to the Terraria source code, here are the parameters for the Eye of Cthulhu in Expert Mode.
    base_health = 3640
    # The health bonus per additional player for the Eye of Cthulhu is a specific 33%.
    player_health_bonus = 0.33

    # The formula for calculating the total health (HP) is:
    # HP = BaseHP * (1 + HealthBonus * (NumberOfPlayers - 1))

    print("The formula for the Eye of Cthulhu's health in Expert Mode multiplayer is based on the number of players.")
    print("\nHere are the components of the equation:")
    print(f"Base Health: {base_health}")
    print(f"Per-Player Health Bonus: {player_health_bonus} (or {player_health_bonus * 100}%)")

    print("\nThe final equation is:")
    print(f"Total Health = {base_health} * (1 + {player_health_bonus} * (NumberOfPlayers - 1))")

    # Let's analyze the formula as the number of players (N) approaches infinity.
    # Health(N) = 3640 * (1 + 0.33*N - 0.33)
    # Health(N) = 3640 * (0.67 + 0.33*N)
    # This simplifies to a linear equation: Health(N) = 2438.8 + 1201.2 * N

    print("\nThis formula is a linear function of the number of players.")
    print("As the number of players increases, the total health increases without any upper limit or cap.")
    print("Therefore, as the number of players approaches infinity, the theoretical health also approaches infinity.")
    print("\nResult: The theoretical maximum health is infinite.")

if __name__ == "__main__":
    calculate_eoc_health()