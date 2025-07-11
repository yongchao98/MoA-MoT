def calculate_terraria_max_boss_health():
    """
    Calculates the theoretical maximum health for the Eye of Cthulhu in Expert Mode.

    This is based on the game's health scaling formula and its internal player cap
    for health calculations.
    """
    # Step 1: Define the base parameters for the Eye of Cthulhu
    base_expert_health = 3640
    per_player_increase_factor = 0.35
    player_cap_for_scaling = 255

    # The problem asks for the theoretical maximum, which is not infinite but
    # capped by the game's code at 255 players for the scaling calculation.

    # Step 2: Calculate the number of "additional" players for the formula
    # The formula is Health = Base * (1 + Factor * (Players - 1))
    additional_players = player_cap_for_scaling - 1

    # Step 3: Calculate the total health multiplier
    health_multiplier = 1 + (per_player_increase_factor * additional_players)

    # Step 4: Calculate the final maximum health
    max_health = base_expert_health * health_multiplier

    # Step 5: Print the equation with all the numbers and the final result.
    # The result is cast to an integer as health points are whole numbers.
    print("Equation for the Eye of Cthulhu's theoretical maximum health in Expert Mode:")
    print(
        f"{base_expert_health} * (1 + {per_player_increase_factor} * ({player_cap_for_scaling} - 1)) = {int(max_health)}"
    )

if __name__ == "__main__":
    calculate_terraria_max_boss_health()

<<<327236>>>