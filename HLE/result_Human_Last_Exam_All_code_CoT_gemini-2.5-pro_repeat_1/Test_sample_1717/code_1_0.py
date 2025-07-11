def calculate_eoc_max_health():
    """
    Calculates the theoretical maximum health of the Eye of Cthulhu in Terraria's Expert Mode.

    In Expert Mode, boss health scales with the number of players. While the prompt asks for
    an infinite number of players, the game's health scaling is capped at the value for 256 players.
    """

    # Base health for the Eye of Cthulhu in Expert Mode.
    base_health = 3640

    # Health bonus per player (after the first) in Expert Mode.
    player_health_bonus = 0.35

    # The game caps the player count for health scaling at 256.
    max_players_for_scaling = 256

    # The health scaling formula is: BaseHealth * (1 + PlayerBonus * (NumberOfPlayers - 1))
    # We substitute the values to find the maximum health.
    
    # Calculate the multiplier
    health_multiplier = 1 + player_health_bonus * (max_players_for_scaling - 1)
    
    # Calculate the final health
    max_health = base_health * health_multiplier

    # Print the explanation and the final equation
    print("The formula for boss health scaling in Expert Mode is: BaseHealth * (1 + PlayerBonus * (NumberOfPlayers - 1))")
    print(f"The calculation is based on the game's maximum player scaling cap of {max_players_for_scaling} players.")
    print("\nFinal Equation:")
    # Using int() for health values as the game uses integers for health.
    print(f"{base_health} * (1 + {player_health_bonus} * ({max_players_for_scaling} - 1)) = {int(max_health)}")

calculate_eoc_max_health()