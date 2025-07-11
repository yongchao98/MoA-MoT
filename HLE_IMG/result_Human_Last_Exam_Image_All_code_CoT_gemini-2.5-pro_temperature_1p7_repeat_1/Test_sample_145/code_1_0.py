def solve_ant_mound_age():
    """
    This script analyzes the provided ecological information to determine the age of two ant mounds.
    """
    age_left_ecosystem = 20  # years since sagebrush was seeded
    age_right_ecosystem = 15 # years since sagebrush was seeded

    print("Analyzing the ant mounds based on the surrounding vegetation.")
    print("-" * 50)

    # --- Analysis of the Right Mound ---
    print("Step 1: Analyze the mound in the right diagram (15-year-old ecosystem).")
    print(f"The sagebrush in this ecosystem was seeded {age_right_ecosystem} years ago.")
    print("Observation: Sagebrush plants are seen growing on the ant mound.")
    print("Reasoning: Harvester ants clear vegetation. If established plants are present on the mound, the ant colony must be younger than the plants, as it has not had enough time to clear them.")
    print(f"Conclusion: The age of the right mound must be less than the age of the sagebrush, so its age is < {age_right_ecosystem} years.")
    print("-" * 50)

    # --- Analysis of the Left Mound ---
    print("Step 2: Analyze the mound in the left diagram (20-year-old ecosystem).")
    print(f"The sagebrush in this ecosystem was seeded {age_left_ecosystem} years ago.")
    print("Observation: The area of the mound is completely clear of sagebrush.")
    print("Reasoning: The absence of vegetation indicates a mature ant colony that has successfully cleared its territory over time.")
    print("Comparing this to the right diagram, this colony is clearly older and more established.")
    print(f"An age in the range of 15 to {age_left_ecosystem} years is consistent with a mature colony in a {age_left_ecosystem}-year-old rehabilitated site.")
    print(f"Conclusion: The age of the left mound is estimated to be between 15 and {age_left_ecosystem} years.")
    print("-" * 50)

    # --- Final Conclusion ---
    print("Step 3: Combine the conclusions.")
    print("Left mound age: 15-20 years")
    print("Right mound age: < 15 years")
    print("This corresponds to option C.")

solve_ant_mound_age()
<<<C>>>