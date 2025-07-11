def find_arabesque_types():
    """
    Analyzes the Vaganova arabesque positions to find which ones have the forward arm
    on the opposite side of the lifted leg.
    """
    
    # In ballet, "opposite" arm/leg is known as "opposition," while "same side" is not.
    # We analyze each arabesque based on this rule. Let's assume the dancer is standing
    # on their right leg, so the left leg is the one lifted behind.
    
    print("To answer the question, we analyze the Vaganova arabesque positions:")
    print("-" * 60)
    
    # First Arabesque
    arabesque_1 = "First"
    lifted_leg_side_1 = "Left"
    forward_arm_side_1 = "Right" # Arm on same side as SUPPORTING leg is forward
    print(f"In {arabesque_1} Arabesque:")
    print(f"  - If the lifted leg is the {lifted_leg_side_1} leg...")
    print(f"  - The forward arm is the {forward_arm_side_1} arm.")
    print(f"  - Conclusion: The forward arm is on the OPPOSITE side of the lifted leg.")
    
    print("-" * 60)

    # Third Arabesque
    arabesque_2 = "Third"
    lifted_leg_side_2 = "Left"
    # In Third Arabesque, both arms are forward. The one corresponding to the supporting leg is lower.
    # This lower arm is the primary focus for the position's "opposition".
    forward_arm_side_2 = "Right" 
    print(f"In {arabesque_2} Arabesque:")
    print(f"  - If the lifted leg is the {lifted_leg_side_2} leg...")
    print(f"  - The arm on the same side as the supporting leg ({forward_arm_side_2}) is forward and lower.")
    print(f"  - Conclusion: This primary forward arm is on the OPPOSITE side of the lifted leg.")
    
    print("-" * 60)

    print("\nBased on this analysis, the two types of arabesque that meet the criteria are:")
    print(f"The {arabesque_1} and {arabesque_2}.")

find_arabesque_types()