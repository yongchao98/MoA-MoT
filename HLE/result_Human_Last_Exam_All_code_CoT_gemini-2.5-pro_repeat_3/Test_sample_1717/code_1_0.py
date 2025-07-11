def solve_terraria_health_limit():
    """
    This script calculates and explains the theoretical maximum health of the
    Eye of Cthulhu in Terraria's Expert Mode with an infinite number of players.
    """
    
    # Define the constants from the game's source code for Expert Mode
    base_expert_health = 3640
    player_scaling_factor = 0.35

    print("Step 1: Identify the base values for the health calculation.")
    print(f"The base health for the Eye of Cthulhu in Expert Mode (1 player) is: {base_expert_health}")
    print(f"The health scaling factor per additional player is: {player_scaling_factor} (or {player_scaling_factor * 100}%)")
    print("-" * 50)

    print("Step 2: State the health scaling formula.")
    print("The formula for total health based on the number of players (n) is:")
    # The final code outputs each number in the final equation, as requested.
    print(f"Total Health = {base_expert_health} * (1 + {player_scaling_factor} * (n - 1))")
    print("-" * 50)
    
    print("Step 3: Simplify the formula to analyze its behavior.")
    print("Let's expand the equation:")
    print(f"Total Health = {base_expert_health} * (1 + {player_scaling_factor}n - {player_scaling_factor})")
    print(f"Total Health = {base_expert_health} * ({(1 - player_scaling_factor):.2f} + {player_scaling_factor}n)")
    
    # Perform the final multiplications
    constant_term = base_expert_health * (1 - player_scaling_factor)
    n_coefficient = base_expert_health * player_scaling_factor
    
    print("Total Health = " + str(base_expert_health) + " * " + f"{(1 - player_scaling_factor):.2f}" + " + " + str(base_expert_health) + " * " + str(player_scaling_factor) + "n")
    print(f"Total Health = {constant_term:.0f} + {n_coefficient:.0f}n")
    print("-" * 50)

    print("Step 4: Draw the final conclusion.")
    print("The health of the Eye of Cthulhu is a linear function of 'n', the number of players.")
    print("As 'n' approaches infinity, the total health also increases without any upper bound.")
    print("\nTherefore, the theoretical maximum health is infinite.")

solve_terraria_health_limit()