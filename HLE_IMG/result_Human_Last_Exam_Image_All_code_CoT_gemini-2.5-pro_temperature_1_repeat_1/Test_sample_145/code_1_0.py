def solve_ant_mound_age():
    """
    This script explains the reasoning to determine the age of the two ant mounds.
    """
    
    # Information about the ecosystems
    ecosystem_left_age = 20  # years
    ecosystem_right_age = 15 # years

    print("Step 1: Analyze the ecosystem on the right.")
    print(f"The sagebrush in the right ecosystem was seeded {ecosystem_right_age} years ago.")
    print("Observation: The image shows sagebrush plants growing within the cleared area of the ant mound.")
    print("Inference: This means the ant colony is younger than the plants and has not yet had enough time to clear them.")
    print("Conclusion: Therefore, the age of the mound on the right is < 15 years.")
    print("-" * 20)

    print("Step 2: Analyze the ecosystem on the left.")
    print(f"The sagebrush in the left ecosystem was seeded {ecosystem_left_age} years ago.")
    print("Observation: The image shows a well-defined cleared disk around the ant mound with no plants inside it.")
    print("Inference: This suggests a mature colony that has been established long enough to clear the surrounding vegetation.")
    print("Comparison: Since a mound < 15 years old (on the right) has NOT cleared the plants, the mound on the left must be older than 15 years to have achieved this.")
    print(f"Conclusion: The mound's age must be less than the ecosystem's age ({ecosystem_left_age} years) but greater than 15 years.")
    print(f"Therefore, the age of the mound on the left is between 15 and 20 years.")
    print("-" * 20)
    
    print("Final Answer Summary:")
    print("Left Mound Age: 15-20 years")
    print("Right Mound Age: < 15 years")
    print("\nThis corresponds to answer choice C.")

solve_ant_mound_age()