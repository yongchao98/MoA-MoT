def solve_ant_mound_age():
    """
    This function analyzes the provided information about two ant mounds
    to determine their approximate ages.
    """

    # Information from the problem description
    age_of_left_ecosystem = 20  # years
    age_of_right_ecosystem = 15 # years

    # Analysis of the Right Mound
    # Observation: The mound is actively clearing sagebrush that is up to 15 years old.
    # Reasoning: For the ants to be clearing existing plants, the ant colony must be
    # younger than the plants.
    age_conclusion_right = f"< {age_of_right_ecosystem}"
    
    # Analysis of the Left Mound
    # Observation: The mound has a large, stable, cleared disc. No plants are being encroached upon.
    # Reasoning: This indicates a mature colony that has been established for a long time,
    # likely for a significant portion of the ecosystem's 20-year history.
    # The mound and plants seem to be in equilibrium.
    age_conclusion_left = f"{age_of_right_ecosystem}-{age_of_left_ecosystem}" # i.e., 15-20

    # Print the step-by-step reasoning
    print("Step 1: Analyze the ecosystem on the right.")
    print(f"The sagebrush in this ecosystem is up to {age_of_right_ecosystem} years old.")
    print("The ant mound is shown clearing these established plants.")
    print("This means the ant colony must be younger than the plants.")
    print(f"Conclusion for Right Mound: Age is {age_conclusion_right} years.\n")

    print("Step 2: Analyze the ecosystem on the left.")
    print(f"The sagebrush in this ecosystem is up to {age_of_left_ecosystem} years old.")
    print("The ant mound has a large, stable cleared area, suggesting it is a mature colony.")
    print("This indicates its age is a significant fraction of the ecosystem's age.")
    print(f"Conclusion for Left Mound: Age is approximately in the {age_conclusion_left} year range.\n")

    print("Step 3: Combine the conclusions.")
    print(f"The age of the left mound is estimated to be {age_conclusion_left} years.")
    print(f"The age of the right mound is estimated to be {age_conclusion_right} years.")
    print("This corresponds to answer choice C.")

solve_ant_mound_age()
<<<C>>>