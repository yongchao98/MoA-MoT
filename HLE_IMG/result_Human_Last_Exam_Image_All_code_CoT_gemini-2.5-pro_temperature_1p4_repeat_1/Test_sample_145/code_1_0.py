def solve_ant_mound_age():
    """
    Determines the age of two ant mounds based on ecological context.
    """
    # Ages of the rehabilitated ecosystems
    ecosystem_age_left = 20
    ecosystem_age_right = 15

    print("Analyzing the ant mound ages based on the provided diagrams and information:")
    print("-" * 70)

    # Reasoning for the right mound
    print(f"1. Right Ecosystem (Seeded {ecosystem_age_right} years ago):")
    print("   - Observation: Sagebrush plants are seen growing within the mound's boundary.")
    print("   - Ecological Principle: Harvester ants (Pogonomyrmex) clear vegetation around their mounds.")
    print(f"   - Deduction: The colony is too young to have cleared the established {ecosystem_age_right}-year-old sagebrush.")
    print(f"   - Conclusion: The age of the right mound is < {ecosystem_age_right} years.")
    print("\n")

    # Reasoning for the left mound
    print(f"2. Left Ecosystem (Seeded {ecosystem_age_left} years ago):")
    print("   - Observation: The mound is a well-defined, clear disk with no vegetation inside.")
    print("   - Deduction: The colony is mature and has had enough time to clear all vegetation in its territory, even plants that could be up to 20 years old.")
    print("   - Comparison: This mound is visibly more established and therefore older than the one on the right.")
    print(f"   - Conclusion: Since the left mound is older than the right mound (<{ecosystem_age_right} years) but younger than its ecosystem ({ecosystem_age_left} years), its age is between {ecosystem_age_right} and {ecosystem_age_left} years.")
    print("\n")

    # Final summary
    print("Summary:")
    print(f"  - Age of Left Mound: 15-20 years")
    print(f"  - Age of Right Mound: <15 years")
    print("\nThis corresponds to answer choice C.")

solve_ant_mound_age()