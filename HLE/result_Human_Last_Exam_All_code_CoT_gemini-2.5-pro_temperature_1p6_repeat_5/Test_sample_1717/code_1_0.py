import math

def calculate_terraria_max_health():
    """
    Calculates and explains the theoretical maximum health of the Eye of Cthulhu
    in Terraria's Expert Mode with an infinite number of players.
    """
    
    # Step 1: Define the base parameters from the game's source code.
    base_health_eoc_expert = 3640
    per_player_multiplier = 0.35  # In Expert Mode, health increases by 35% for each player after the first.
    
    # The health formula is: Health = BaseHealth * (1 + Multiplier * (NumberOfPlayers - 1))
    
    print("Terraria Expert Mode Health Calculation for Eye of Cthulhu\n")
    print(f"Base Health: {base_health_eoc_expert}")
    print(f"Per-Player Multiplier: {per_player_multiplier*100}%")
    print("Formula: Health = BaseHealth * (1 + 0.35 * (N - 1))")
    print("-" * 50)
    
    # Step 2: Analyze the "infinite players" constraint.
    print("With an infinite number of players (N), the health calculated by the formula would also be infinite.")
    print("However, Terraria's source code reveals a technical limitation.\n")
    
    # Step 3: Introduce the engine's hard cap.
    # NPC health is stored as a 32-bit signed integer (int), which has a maximum value.
    health_cap = 2147483647
    
    print(f"NPC health is stored in a 32-bit signed integer, which has a maximum value of {health_cap:,}.")
    print("Any calculated health value exceeding this cap is clamped to this maximum value.")
    print("Therefore, the theoretical maximum health is this hard-coded limit.\n")
    
    # Step 4: Show the final conclusion and the numbers involved.
    # We can show the equation with the capped value as the result.
    # MaxHealth = BaseHealth * (1 + 0.35 * (N - 1))
    
    print("Final Equation:")
    print("The final health is the lesser of the two values: the result of the formula or the integer limit.")
    print("Final Health = min(FormulaResult, EngineCap)")
    print(f"So, the final theoretical maximum health is the engine cap itself.")
    
    # The prompt asks to output each number in the final equation.
    # We can present it as: "EngineCap = The final result of the equation under max conditions"
    # To satisfy this, let's find the number of players (N) where the cap is reached.
    # health_cap = base_health_eoc_expert * (1 + per_player_multiplier * (N-1))
    # (health_cap / base_health_eoc_expert) = 1 + per_player_multiplier * (N-1)
    # ((health_cap / base_health_eoc_expert) - 1) / per_player_multiplier = N - 1
    # N = (((health_cap / base_health_eoc_expert) - 1) / per_player_multiplier) + 1
    
    players_to_reach_cap = (((health_cap / base_health_eoc_expert) - 1) / per_player_multiplier) + 1
    
    print("\nTo demonstrate, the number of players required to hit this cap is roughly:")
    print(f"N = (( {health_cap} / {base_health_eoc_expert} ) - 1) / {per_player_multiplier} + 1 ≈ {math.ceil(players_to_reach_cap):,}")

    print("\n--- FINAL ANSWER ---")
    final_health = health_cap
    print(f"The theoretical maximum health of the Eye of Cthulhu is {final_health:,}.")


calculate_terraria_max_health()
<<<2147483647>>>