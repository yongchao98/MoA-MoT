import math

def calculate_max_health():
    """
    Calculates the theoretical maximum health of the Eye of Cthulhu in Expert Mode.
    
    In Terraria's source code, the health for NPCs is stored as a 32-bit signed integer.
    While the health scaling formula with infinite players would mathematically result in infinite health,
    the actual health is capped by the maximum value this data type can hold.
    """
    
    # Step 1: Define the base parameters for the Eye of Cthulhu in Expert Mode.
    base_health_expert = 3640
    health_increase_per_player = 0.35 # 35%
    
    print(f"The Eye of Cthulhu's base health in Expert Mode is {base_health_expert}.")
    print(f"Health increases by {health_increase_per_player * 100}% for each player after the first.")
    print("\nThe health formula is: Health = BaseHealth * (1 + 0.35 * (NumberOfPlayers - 1))")
    print("This formula increases linearly and has no mathematical maximum.")
    
    # Step 2: Consider the technical limitation from the game's source code.
    # NPC health is stored in a 32-bit signed integer (int in C#).
    # The maximum value for a 32-bit signed integer is 2^31 - 1.
    
    base = 2
    exponent = 31
    subtrahend = 1
    
    max_health_value = (base ** exponent) - subtrahend
    
    print("\nHowever, the game's code imposes a technical limit due to the 32-bit signed integer data type used for health.")
    print(f"The final equation for the maximum possible health is derived from this limit: {base} ^ {exponent} - {subtrahend}")
    print(f"\nThe theoretical maximum health is: {max_health_value:,}")

calculate_max_health()