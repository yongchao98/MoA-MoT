def explain_mound_ages():
    """
    This script explains the reasoning for determining the ant mound ages
    based on the provided ecological diagrams and information.
    """

    left_ecosystem_age = 20
    right_ecosystem_age = 15
    
    print("### Analysis of the Ecosystems ###")
    
    print("\n--- Left Ecosystem (Seeded {} years ago) ---".format(left_ecosystem_age))
    print("Observation: The ant mound is situated between established sagebrush plants.")
    print("Inference: The plants were present before the ant mound was created.")
    print("Conclusion: The age of the left mound is LESS THAN the age of the ecosystem.")
    print("=> Left Mound Age < {} years".format(left_ecosystem_age))
    
    print("\n--- Right Ecosystem (Seeded {} years ago) ---".format(right_ecosystem_age))
    print("Observation: The sagebrush plants grow around the mound's clearing, which interrupts the plant grid.")
    print("Inference: The mound's clearing existed before the area was seeded.")
    print("Conclusion: The age of the right mound is GREATER THAN the age of the seeding.")
    print("=> Right Mound Age > {} years".format(right_ecosystem_age))
    
    print("\n### Final Determination ###")
    print("Based on the analysis, the Left Mound is < {} years old and the Right Mound is > {} years old.".format(left_ecosystem_age, right_ecosystem_age))
    print("To refine this, we can assume a regional event ~20 years ago prompted both rehabilitations.")
    print("This puts an upper limit on the right mound's age, making it also < 20 years.")
    print("Therefore, the age range for the Right Mound is between {} and {} years.".format(right_ecosystem_age, left_ecosystem_age))
    print("The choices suggest the Left Mound's age is < 15 years, which fits our < 20 years conclusion.")
    
    print("\nFinal Answer Logic:")
    print("Left Mound Age: < 15 years")
    print("Right Mound Age: 15-20 years")
    print("This corresponds to answer choice D.")

explain_mound_ages()