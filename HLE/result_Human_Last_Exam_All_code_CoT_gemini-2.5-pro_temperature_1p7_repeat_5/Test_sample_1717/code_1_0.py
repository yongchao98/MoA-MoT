def calculate_eoc_health():
    """
    Calculates and explains the theoretical maximum health of the Eye of Cthulhu
    in Terraria Expert Mode with an infinite number of players.
    """
    
    # Step 1: Define the base health for the Eye of Cthulhu in Expert Mode.
    base_health_expert = 3640
    
    # Step 2: Define the health scaling factor per player in Expert Mode.
    # Health increases by 35% of the base health for each player after the first.
    scaling_per_player = 0.35
    
    # Calculate the flat amount of health added per additional player.
    health_increase_per_player = int(base_health_expert * scaling_per_player)
    
    print("Determining the theoretical maximum health for the Eye of Cthulhu in Expert Mode.")
    print("-" * 70)
    
    # Step 3: Present the formula for N players.
    print(f"The base health for one player is: {base_health_expert}")
    print(f"For each additional player, health increases by {scaling_per_player:.0%}, which is a flat value of {health_increase_per_player} HP.")
    print("\nThis gives us the general health equation, where 'N' is the number of players:")
    
    # The prompt requires outputting each number in the final equation.
    print("\nFinal Equation:")
    print(f"Total Health = {base_health_expert} + ({health_increase_per_player} * (N - 1))")
    
    # Step 4: Analyze the formula for an infinite number of players.
    print("\nAnalysis for an infinite number of players:")
    print("The health of the Eye of Cthulhu increases linearly with the number of players ('N').")
    print("As the value of 'N' grows larger, the Total Health continues to increase without any upper bound.")
    print("\nTherefore, as the number of players approaches infinity, the boss's health also approaches infinity.")
    print("There is no finite theoretical maximum health value.")

calculate_eoc_health()