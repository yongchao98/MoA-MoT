def calculate_eoc_health():
    """
    Calculates and explains the theoretical maximum health of the Eye of Cthulhu
    in Terraria's Expert Mode with an infinite number of players.
    """

    # Step 1: Define the base constants from the game's code.
    base_health = 2800
    expert_mode_health_solo = 3640

    # The Expert Mode multiplier is the ratio of solo expert health to base health.
    expert_multiplier = expert_mode_health_solo / base_health

    # Step 2: Define the multiplayer health scaling factor per player (after the first).
    player_scaling_factor = 0.35

    print("The formula for boss health in Expert Mode multiplayer is:")
    print("Total Health = Base Health * Expert Multiplier * (1 + Player Scaling Factor * (n - 1))\n")

    print("For the Eye of Cthulhu:")
    print(f"Base Health = {base_health}")
    print(f"Expert Multiplier = {expert_mode_health_solo} / {base_health} = {expert_multiplier}")
    print(f"Player Scaling Factor = {player_scaling_factor}\n")

    print("Let 'n' be the number of players.")
    print("The formula is: Health(n) = 2800 * 1.3 * (1 + 0.35 * (n - 1))")
    print("Simplifying the expression:")
    print("Health(n) = 3640 * (1 + 0.35n - 0.35)")
    print("Health(n) = 3640 * (0.65 + 0.35n)")
    print("Health(n) = (3640 * 0.65) + (3640 * 0.35)n\n")

    # Step 3: Calculate the final simplified linear equation terms.
    # This is the y = mx + c form, where x is n.
    health_intercept = base_health * expert_multiplier * (1 - player_scaling_factor)
    health_slope = base_health * expert_multiplier * player_scaling_factor

    print("This gives us the final simplified linear equation:")
    print(f"Health(n) = {int(health_slope)}n + {int(health_intercept)}\n")

    print("As you can see, the health of the Eye of Cthulhu is a linear function of 'n', the number of players.")
    print("As 'n' approaches infinity, the total health also approaches infinity.")
    print("Therefore, there is no finite theoretical maximum health.")

calculate_eoc_health()