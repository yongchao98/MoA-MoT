def solve_ant_mound_age():
    """
    This function explains the reasoning to determine the age of the ant mounds.
    """
    
    # Information about the ecosystems
    left_ecosystem_age = 20  # years since seeding
    right_ecosystem_age = 15 # years since seeding

    print("Step 1: Analyze the ecosystem on the right.")
    print(f"The sagebrush was seeded {right_ecosystem_age} years ago.")
    print("The diagram shows that sagebrush plants are still growing on the edge of the ant mound's clearing.")
    print("This indicates the ant colony is not yet mature enough to have fully cleared the established vegetation.")
    print(f"Therefore, the age of the right ant mound must be less than {right_ecosystem_age} years.")
    print("\nConclusion for Right Mound: Age < 15 years.\n")

    print("Step 2: Analyze the ecosystem on the left.")
    print(f"The sagebrush was seeded {left_ecosystem_age} years ago.")
    print("The diagram shows a complete, well-defined clearing with no sagebrush growing near the mound.")
    print("This implies the ant colony is mature and has been established long enough to clear its territory.")
    print(f"In a {left_ecosystem_age}-year-old ecosystem, a mature colony that has created such a clearing would likely be in the 15 to 20 year age range.")
    print("\nConclusion for Left Mound: Age is approximately 15-20 years.\n")

    print("Step 3: Combine the conclusions.")
    print("Left Mound Age: 15-20 years")
    print("Right Mound Age: < 15 years")
    print("This corresponds to answer choice C.")

solve_ant_mound_age()