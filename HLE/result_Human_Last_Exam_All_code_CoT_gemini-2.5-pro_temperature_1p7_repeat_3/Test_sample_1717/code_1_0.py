def solve_terraria_health_cap():
    """
    Calculates and explains the theoretical maximum health for the Eye of Cthulhu
    in Terraria's Expert Mode with an infinite number of players.
    """
    
    # Step 1: Define the base parameters from the game's code.
    base_health = 2800
    expert_mode_multiplier = 2.0
    multiplayer_health_increase_per_player = 0.35 # This is a 35% increase per player after the first one.

    # Step 2: Calculate the base health in Expert Mode.
    expert_base_health = base_health * expert_mode_multiplier
    
    # Step 3: Explain the formula and the process.
    print("To find the theoretical maximum health of the Eye of Cthulhu, we follow these steps:")
    print("---------------------------------------------------------------------------------")
    print(f"1. The base health of the Eye of Cthulhu is {base_health}.")
    print(f"2. In Expert Mode, its health is multiplied by {expert_mode_multiplier}.")
    print(f"   The resulting single-player Expert Mode health is: {base_health} * {expert_mode_multiplier} = {int(expert_base_health)}")
    print("\n3. In multiplayer, health increases for each player after the first.")
    print("   The formula is: Total Health = Expert_Base_Health * (1 + Scaling_Factor * (Number_Of_Players - 1))")
    print(f"   Substituting the values: Total Health = {int(expert_base_health)} * (1 + {multiplayer_health_increase_per_player} * (P - 1)), where P is the number of players.")

    # Step 4: Analyze the limit and technical constraints.
    print("\n4. With an infinite number of players (as P -> infinity), the health calculated by the formula also approaches infinity.")
    print("   However, games have technical limitations based on their code.")
    print("\n5. In Terraria's source code, an NPC's health is stored in a 32-bit signed integer.")
    
    # The maximum value for a 32-bit signed integer.
    max_int_32_value = 2147483647
    
    print(f"   A 32-bit signed integer has a maximum possible value of {max_int_32_value}.")
    print("   Therefore, the game engine cannot store a health value higher than this, regardless of the formula's output.")
    
    # Final Conclusion
    print("\n---------------------------------------------------------------------------------")
    print("The theoretical maximum health is not infinite, but capped by this data type limit.")
    print(f"Final Answer: The theoretical maximum health is {max_int_32_value}.")


solve_terraria_health_cap()