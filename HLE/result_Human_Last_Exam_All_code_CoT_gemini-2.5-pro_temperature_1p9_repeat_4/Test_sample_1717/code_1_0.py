import math

def calculate_eoc_health():
    """
    Calculates and explains the theoretical maximum health of the Eye of Cthulhu
    in Terraria Expert Mode with an infinite number of players.
    """
    # Step 1: Define the base constants from the game's code for Expert Mode.
    # Base health for the Eye of Cthulhu in Expert Mode.
    base_health = 3640
    # Health increase per additional player in Expert Mode (35%).
    player_increase_percentage = 0.35

    print("--- Terraria Health Calculation for Eye of Cthulhu (Expert Mode) ---")

    # Step 2: State the general formula for boss health scaling in multiplayer.
    print("\nThe general formula for a boss's health in Expert Mode multiplayer is:")
    print("Total Health = Base Health * (1 + (Number of Players - 1) * Player Increase)")
    print("\nLet 'N' be the number of players.")
    print(f"The formula for the Eye of Cthulhu is:")
    print(f"Total Health = {base_health} * (1 + (N - 1) * {player_increase_percentage})")

    # Step 3: Simplify the equation step-by-step.
    print("\nSimplifying the equation:")
    print(f"1. Distribute the percentage: Total Health = {base_health} * (1 + {player_increase_percentage}*N - {player_increase_percentage})")
    
    one_minus_increase = 1 - player_increase_percentage
    print(f"2. Combine constants inside parentheses: Total Health = {base_health} * ({one_minus_increase:.2f} + {player_increase_percentage}*N)")
    
    # Calculate the final constant and multiplier
    constant_term = base_health * one_minus_increase
    multiplier_term = base_health * player_increase_percentage
    
    print(f"3. Distribute the base health: Total Health = ({base_health} * {one_minus_increase:.2f}) + ({base_health} * {player_increase_percentage} * N)")

    # Step 4: Present the final simplified linear equation.
    print("\nThe final simplified equation for the Eye of Cthulhu's health is:")
    print(f"Total Health = {math.ceil(constant_term)} + {math.ceil(multiplier_term)} * N")

    # Step 5: Conclude based on the nature of the equation for an infinite number of players.
    print("\nThis equation shows that the Total Health increases linearly with 'N', the number of players.")
    print("As 'N' approaches infinity, the Total Health will also increase without any upper limit.")
    print("\nTherefore, the theoretical maximum health is infinite.")


if __name__ == '__main__':
    calculate_eoc_health()