import math

def calculate_max_health():
    """
    Calculates the theoretical maximum health of the Eye of Cthulhu in Terraria's Expert Mode.
    The calculation is based on an infinite number of players, which is limited by the
    game's data type for storing health.
    """

    # Step 1: Define the base health values from the game's data.
    # We will use the health of the Eye of Cthulhu's first phase for calculations.
    classic_base_health = 2800
    expert_single_player_health = 3640

    # Step 2: Define the health scaling mechanic for additional players in Expert Mode.
    # Each player after the first adds 50% of the boss's Classic Mode health.
    per_player_scaling_factor = 0.50
    health_per_additional_player = int(classic_base_health * per_player_scaling_factor)

    # Step 3: Explain the health formula for N players.
    print("--- Terraria Boss Health Calculation ---")
    print("\nStep 1: The health scaling formula for a boss in Expert Mode multiplayer is:")
    print("Total Health = (Single-Player Expert Health) + (Health Per Additional Player) * (Number of Players - 1)\n")

    # Step 4: Show the specific values for the Eye of Cthulhu in the equation.
    print("Step 2: Plugging in the values for the Eye of Cthulhu gives the equation:")
    print(f"Total Health = {expert_single_player_health} + ({classic_base_health} * {per_player_scaling_factor}) * (N - 1)")
    print(f"Total Health = {expert_single_player_health} + {health_per_additional_player} * (N - 1)\n")
    
    # Step 5: Analyze the concept of "infinite players".
    print("Step 3: Analyzing the 'infinite players' condition.")
    print("With the formula above, as the number of players (N) approaches infinity, the total health also approaches infinity.")
    print("However, Terraria's source code stores NPC health as a 32-bit signed integer.\n")

    # Step 6: Identify the technical limit.
    # The maximum value for a 32-bit signed integer is 2^31 - 1.
    max_int_32_bit = 2**31 - 1

    print("Step 4: The technical limitation imposes a maximum health cap.")
    print(f"The maximum value for a 32-bit signed integer is {max_int_32_bit:,}.")
    print("This is the theoretical maximum health for any single NPC in the game.")

calculate_max_health()