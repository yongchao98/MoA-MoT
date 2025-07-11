def solve_ant_mound_age():
    """
    Analyzes the provided information about ant mounds and sagebrush to determine the age of each mound.
    """
    # Ages of the rehabilitated ecosystems
    age_eco_left = 20  # years
    age_eco_right = 15 # years

    print("Analyzing the age of the ant mounds based on ecological principles.")
    print("-" * 30)

    # Analysis of the right mound
    print(f"1. The ecosystem on the right was seeded with sagebrush {age_eco_right} years ago.")
    print("2. The diagram shows sagebrush plants growing *inside* the boundary of the ant mound.")
    print("3. This indicates the sagebrush was established *before* the ant colony grew large enough to clear it.")
    print(f"4. Therefore, the ant mound on the right must be younger than the plants.")
    age_mound_right = f"< {age_eco_right}"
    print(f"Conclusion for the right mound: Age is {age_mound_right} years.")
    print("-" * 30)

    # Analysis of the left mound
    print(f"1. The ecosystem on the left was seeded with sagebrush {age_eco_left} years ago.")
    print("2. The diagram shows a mound with a completely cleared disk, with no sagebrush growing inside.")
    print("3. This indicates a mature ant colony that has prevented sagebrush from establishing in its territory over a long period.")
    print(f"4. The colony must be well-established, likely existing for a significant portion of the ecosystem's {age_eco_left}-year history.")
    age_mound_left = f"{age_eco_right}-{age_eco_left}" # "15-20"
    print(f"Conclusion for the left mound: Age is in the {age_mound_left} year range.")
    print("-" * 30)
    
    # Final Answer
    print(f"Final Answer Summary:")
    print(f"Left Mound Age: {age_mound_left} years")
    print(f"Right Mound Age: {age_mound_right} years")
    print("\nThis corresponds to answer choice C.")

solve_ant_mound_age()