def solve_ant_mound_age():
    """
    Determines the age of two ant mounds based on ecological information.
    """
    # Information provided in the problem
    age_ecosystem_left = 20  # years
    age_ecosystem_right = 15 # years

    # Observations from the diagrams
    # Right diagram: Sagebrush plants are inside the mound's cleared area.
    # Left diagram: The mound's cleared area is free of sagebrush plants.

    print("Analyzing the ant mound in the right ecosystem:")
    print(f"The ecosystem was seeded {age_ecosystem_right} years ago.")
    print("Observation: Sagebrush plants are still present inside the mound's clearing.")
    print("Reasoning: This implies the ant colony is younger than the established plants and has not had enough time to clear them.")
    print(f"Conclusion: The age of the right mound is < {age_ecosystem_right} years.")
    print("-" * 30)

    print("Analyzing the ant mound in the left ecosystem:")
    print(f"The ecosystem was seeded {age_ecosystem_left} years ago.")
    print("Observation: The mound's clearing is complete, with no sagebrush inside.")
    print("Reasoning: This indicates the ant colony is older and more mature, having had enough time to clear all vegetation.")
    print(f"Comparing the two, the left mound must be older than the right mound (<{age_ecosystem_right} years) but cannot be older than the ecosystem itself ({age_ecosystem_left} years).")
    print(f"Conclusion: The age of the left mound is between {age_ecosystem_right} and {age_ecosystem_left} years.")
    print("-" * 30)

    print("Final Answer Summary:")
    print(f"Left Mound Age: {age_ecosystem_right}-{age_ecosystem_left} years")
    print(f"Right Mound Age: <{age_ecosystem_right} years")
    print("This corresponds to option C.")

solve_ant_mound_age()
<<<C>>>