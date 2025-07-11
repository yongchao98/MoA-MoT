import math

def solve_terraria_health_cap():
    """
    Calculates the theoretical maximum health for the Eye of Cthulhu
    in Expert Mode with an infinite number of players, based on Terraria's source code.
    """
    # Step 1: Define the base health from the game's code.
    base_health = 2800

    # Step 2: Define the Expert Mode health multiplier.
    expert_mode_multiplier = 2.0

    # Step 3: Define the health scaling formula for 'P' players.
    # The formula is: Health(P) = (base_health * expert_mode_multiplier) * (1 + 0.35 * P)
    # As the number of players 'P' approaches infinity, the result of this formula also approaches infinity.

    # Step 4: Identify the hard-coded engine limit from the source code.
    # NPC health is capped at the maximum value for a 32-bit signed integer.
    health_cap = int(math.pow(2, 31) - 1)

    # Step 5: The final answer is this cap, as the unbounded calculation will always exceed it.
    # We will print the numbers that lead to this conclusion.
    print("The unbounded health formula is: (Base Health * Expert Multiplier) * (1 + 0.35 * Player_Count)")
    print(f"This would be: ({base_health} * {expert_mode_multiplier}) * (1 + 0.35 * infinity), which tends to infinity.")
    print("\nHowever, the game's source code caps all NPC health at the maximum value for a 32-bit signed integer.")
    print("\nThe final equation is the value of this engine cap:")
    
    # Final Output
    print(f"Theoretical Maximum Health = {health_cap}")

solve_terraria_health_cap()