def calculate_terraria_boss_health():
    """
    Calculates the theoretical maximum health of the Eye of Cthulhu in Terraria's Expert Mode.

    According to the Terraria source code, boss health in Expert Mode scales with the number of active players.
    While the user asks for an "infinite" number of players, the game's code has a practical limit based on
    its player array size, which is 255. This script calculates the health based on that limit.
    """

    # 1. Base health of Eye of Cthulhu in Expert Mode.
    base_health = 3640

    # 2. The maximum number of players the game's player array can hold.
    max_players = 255

    # 3. The health scaling formula is: BaseHealth * (1 + 0.35 * (NumberOfPlayers - 1))
    # We calculate each part of the equation step-by-step.

    players_minus_one = max_players - 1
    scaling_per_player = 0.35
    total_scaling_factor = scaling_per_player * players_minus_one
    final_multiplier = 1 + total_scaling_factor
    
    # The game calculates health using floating-point numbers and then casts to an integer.
    max_health = int(base_health * final_multiplier)

    # 4. Print the full equation and the final answer.
    print("Health Scaling Formula: Base Health * (1 + 0.35 * (Number of Players - 1))")
    print(f"Based on the game's code, the maximum number of players is {max_players}.")
    print("\nCalculation:")
    print(f"{base_health} * (1 + {scaling_per_player} * ({max_players} - 1))")
    print(f"= {base_health} * (1 + {scaling_per_player} * {players_minus_one})")
    print(f"= {base_health} * (1 + {total_scaling_factor})")
    print(f"= {base_health} * {final_multiplier}")
    print(f"= {max_health}")


if __name__ == "__main__":
    calculate_terraria_boss_health()