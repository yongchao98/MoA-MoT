def calculate_terraria_boss_health():
    """
    Calculates the theoretical maximum health of the Eye of Cthulhu in Terraria's Expert Mode.
    """
    # Step 1: Define the constants based on game mechanics.
    base_health = 3640  # Eye of Cthulhu's base health in Expert Mode.
    health_increase_per_player = 0.35  # 35% health increase per additional player.
    max_players_for_scaling = 255  # The game caps health scaling at 255 players.

    # Step 2: The number of players that contribute to the health increase
    # is the total number of players minus the first one.
    extra_players = max_players_for_scaling - 1

    # Step 3: Calculate the total health multiplier.
    total_multiplier = 1 + (health_increase_per_player * extra_players)

    # Step 4: Calculate the final maximum health.
    max_health = base_health * total_multiplier

    # Print the explanation and the final equation with all numbers.
    print("The formula for boss health scaling in Expert Mode is:")
    print("Max Health = Base Health * (1 + Health Increase * (Max Players - 1))")
    print("\nApplying the values for the Eye of Cthulhu:")
    print(f"Max Health = {base_health} * (1 + {health_increase_per_player} * ({max_players_for_scaling} - 1))")
    print(f"Max Health = {base_health} * ({1 + health_increase_per_player * extra_players})")
    print(f"Max Health = {base_health} * {total_multiplier}")
    print(f"\nThe theoretical maximum health of the Eye of Cthulhu is: {int(max_health)}")

calculate_terraria_boss_health()
<<<327236>>>