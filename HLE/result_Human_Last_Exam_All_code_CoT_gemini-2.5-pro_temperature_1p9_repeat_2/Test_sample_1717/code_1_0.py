import math

def calculate_max_health():
    """
    Calculates the theoretical maximum health for an NPC in Terraria.
    """
    # Step 1: Define the base health and scaling factor for context.
    base_health = 3640
    scaling_per_player = 0.35
    health_increase = base_health * scaling_per_player

    print(f"The health scaling formula for the Eye of Cthulhu in Expert Mode is:")
    print(f"Total Health = {base_health} (Base) + {int(health_increase)} * (Number of Players - 1)\n")
    
    print("With an infinite number of players, this formula would approach infinity.")
    print("However, Terraria's engine has a hard limit on NPC health based on its data type.\n")

    # Step 2: Define the components of the actual technical limit.
    # Health is stored as a 32-bit signed integer, which has a maximum value of 2^31 - 1.
    base = 2
    exponent = 31
    subtract = 1
    
    # Step 3: Calculate the maximum value.
    max_health = int(math.pow(base, exponent)) - subtract

    # Step 4: Print the final equation and the result.
    print("The theoretical maximum health is capped by the maximum value of a 32-bit signed integer.")
    print("The equation for this limit is:")
    print(f"{base} ** {exponent} - {subtract} = {max_health}")
    print("\nTherefore, the theoretical maximum health of the Eye of Cthulhu is this value.")

calculate_max_health()