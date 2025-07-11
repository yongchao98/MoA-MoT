def calculate_terraria_boss_health():
    """
    Calculates the theoretical maximum health of the Eye of Cthulhu in Expert Mode.

    According to the Terraria source code, boss health scales with the number of players.
    While the question asks for infinite players, the game's architecture imposes a hard
    limit of 256 players in the Main.player array. This script calculates the health
    based on this technical limitation.
    """

    # Base health for Eye of Cthulhu in Expert Mode
    base_health = 3640

    # The game's engine has a maximum player limit of 256.
    num_players = 256

    # In Expert Mode, boss health increases by 35% for each player after the first.
    health_increase_per_player = 0.35

    # The formula is: BaseHealth * (1 + (NumberOfPlayers - 1) * HealthIncrease)
    multiplier = 1 + (num_players - 1) * health_increase_per_player

    # Final health is always an integer in the game.
    final_health = int(base_health * multiplier)

    print("Terraria Boss Health Calculation (Expert Mode)")
    print("Boss: Eye of Cthulhu")
    print("-" * 40)
    print(f"Base Health: {base_health}")
    print(f"Player Count (Theoretical Max): {num_players}")
    print(f"Health Increase Per Additional Player: {health_increase_per_player:.2f}")
    print("-" * 40)
    print("Final Equation:")
    print(f"{base_health} * (1 + ({num_players} - 1) * {health_increase_per_player}) = {final_health}")


if __name__ == "__main__":
    calculate_terraria_boss_health()