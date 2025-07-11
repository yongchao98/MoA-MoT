def calculate_terraria_boss_health():
    """
    Calculates the theoretical maximum health of the Eye of Cthulhu in Expert Mode.

    According to the Terraria source code:
    1. The base health for the Eye of Cthulhu in Expert Mode is 3640.
    2. Boss health scales with the number of players. The formula is:
       Total Health = Base Health * (1 + 0.35 * (Number of Players - 1))
    3. The "theoretical maximum" based on the source code is limited by the maximum
       number of players the game engine supports, which is 255.
    """
    
    # Base health of Eye of Cthulhu in Expert Mode.
    base_health = 3640
    
    # The maximum number of players supported by the Terraria engine.
    num_players = 255
    
    # The health scaling factor per additional player in Expert Mode.
    scaling_factor = 0.35
    
    # The number of players beyond the first who contribute to health scaling.
    additional_players = num_players - 1
    
    # Calculate the total health using the game's formula.
    total_health = base_health * (1 + scaling_factor * additional_players)

    # Print the equation with all the numbers filled in.
    print(f"Base Health * (1 + Scaling Factor * (Number of Players - 1))")
    print(f"{base_health} * (1 + {scaling_factor} * ({num_players} - 1)) = {int(total_health)}")

calculate_terraria_boss_health()