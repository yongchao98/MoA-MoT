def solve_mound_age():
    """
    This function analyzes the provided ecological information to determine the age of two ant mounds.
    """

    # Information for the ecosystem on the right
    right_ecosystem_age = 15  # years
    right_observation = "Sagebrush plants are growing within the ant mound's boundary."

    # Information for the ecosystem on the left
    left_ecosystem_age = 20  # years
    left_observation = "There is a complete clearing around the ant mound with no sagebrush inside."

    print("Step 1: Analyze the mound on the right.")
    print(f"The ecosystem was seeded {right_ecosystem_age} years ago.")
    print(f"Observation: {right_observation}")
    print("This indicates the ant colony is not yet established enough to clear the existing vegetation.")
    print(f"Conclusion: The right mound must be younger than the plants. Its age is < {right_ecosystem_age} years.\n")

    print("Step 2: Analyze the mound on the left.")
    print(f"The ecosystem was seeded {left_ecosystem_age} years ago.")
    print(f"Observation: {left_observation}")
    print("This shows a mature colony that has successfully cleared its territory.")
    print("Compared to the right mound, this colony must be older and more established.")
    print("An age range of 15-20 years represents a mature colony that could create such a clearing within a 20-year-old ecosystem.\n")

    print("Step 3: Combine the conclusions.")
    print("Age of the left mound: 15-20 years.")
    print("Age of the right mound: < 15 years.")
    print("This corresponds to option C.\n")

solve_mound_age()