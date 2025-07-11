def solve_ecology_problem():
    """
    This function explains the reasoning to determine the age of the ant mounds.
    """
    
    # Information about the ecosystems
    age_of_left_ecosystem = 20  # years
    age_of_right_ecosystem = 15 # years
    
    # --- Analysis of the Right Mound ---
    print("Step 1: Analyze the ecosystem on the right.")
    print(f"The sagebrush in this ecosystem was seeded {age_of_right_ecosystem} years ago.")
    print("Observation: Sagebrush plants are growing on the ant mound.")
    print("Inference: Pogonomyrmex ants clear vegetation. Since plants are present on the mound, the ant colony must be younger than the plants.")
    print(f"Conclusion for Right Mound: The age of the mound is < {age_of_right_ecosystem} years.")
    print("-" * 20)
    
    # --- Analysis of the Left Mound ---
    print("Step 2: Analyze the ecosystem on the left.")
    print(f"The sagebrush in this ecosystem was seeded {age_of_left_ecosystem} years ago.")
    print("Observation: The ant mound is completely clear of any vegetation.")
    print("Inference: The ant colony is mature and established enough to have cleared the area and prevented the sagebrush from growing on it.")
    print(f"This implies the colony has been active for a long time, likely establishing early in the {age_of_left_ecosystem}-year history of the site.")
    print(f"Conclusion for Left Mound: A reasonable age estimate is between 15 and {age_of_left_ecosystem} years.")
    print("-" * 20)

    # --- Final Conclusion ---
    print("Step 3: Combine the conclusions.")
    print("Left mound age: 15-20 years")
    print("Right mound age: <15 years")
    print("This corresponds to answer choice C.")

solve_ecology_problem()