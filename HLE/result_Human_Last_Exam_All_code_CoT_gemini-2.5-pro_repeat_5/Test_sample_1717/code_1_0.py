def calculate_eoc_max_health():
    """
    Calculates and explains the theoretical maximum health of the Eye of Cthulhu
    in Terraria's Expert Mode with an infinite number of players.
    """

    # According to Terraria's source code, the base health for the
    # Eye of Cthulhu in Expert Mode (for one player) is 3640.
    base_health = 3640

    # In Expert Mode, boss health is increased by a certain percentage
    # for each player after the first. For the Eye of Cthulhu, this value is 35%.
    health_scaling_factor = 0.35

    # The formula to calculate health (H) based on the number of players (n) is:
    # H(n) = BaseHealth * (1 + HealthScalingFactor * (n - 1))

    print("The formula for the Eye of Cthulhu's health in Expert Mode is:")
    print("Health = BaseHealth * (1 + ScalingFactor * (NumberOfPlayers - 1))\n")

    print("Based on the game's data, the specific values are:")
    print(f"BaseHealth = {base_health}")
    print(f"ScalingFactor = {health_scaling_factor}\n")

    print("So, the final equation for health (H) with 'n' players is:")
    print(f"H(n) = {base_health} * (1 + {health_scaling_factor} * (n - 1))\n")

    print("This equation simplifies to a linear function: H(n) = 1274 * n + 2366.")
    print("As the number of players 'n' approaches infinity, the result of this linear function also approaches infinity.")
    print("Therefore, there is no finite maximum health for the Eye of Cthulhu in this scenario.")
    print("\nThe theoretical maximum health is infinite.")

if __name__ == "__main__":
    calculate_eoc_max_health()