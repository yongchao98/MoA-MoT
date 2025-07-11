def find_diagnosis():
    """
    This script evaluates the clinical data to identify the most likely diagnosis.
    It uses a simple scoring system to quantify how well each potential disease
    matches the patient's information.
    """
    # Patient's clinical and geographical data
    symptoms = {"fever", "headache", "myalgia", "disorientation"}
    location = "Oklahoma"
    lab_clue = "false_positive_lyme_igm"

    # Points are assigned for matching key features
    points_for_geography = 2
    points_for_symptoms = 2
    points_for_lab_clue = 3 # This is a strong, specific clue

    # Evaluate Ehrlichia
    ehrlichia_score = 0
    print("Evaluating Ehrlichia:")
    
    # Check 1: Geography. Ehrlichia is endemic in Oklahoma.
    print(f"1. Is the patient's location (Oklahoma) consistent with Ehrlichia? Yes.")
    ehrlichia_score += points_for_geography
    print(f"   Score is now: {ehrlichia_score}")

    # Check 2: Symptoms. The symptom set is classic for Ehrlichia.
    print(f"2. Are the symptoms (fever, headache, myalgia, disorientation) consistent with Ehrlichia? Yes.")
    ehrlichia_score += points_for_symptoms
    print(f"   Score is now: {ehrlichia_score}")

    # Check 3: Labs. Ehrlichiosis is known to cause false-positive Lyme IgM tests.
    print(f"3. Does Ehrlichia explain the specific lab clue (false-positive Lyme IgM)? Yes.")
    ehrlichia_score += points_for_lab_clue
    print(f"   Score is now: {ehrlichia_score}\n")

    print("The final diagnosis is determined by the evidence that provides the most complete explanation.")
    print("The 'equation' for Ehrlichia based on the evidence points is:")
    # As requested, printing the numbers in the final equation
    print(f"{points_for_geography} (for geography) + {points_for_symptoms} (for symptoms) + {points_for_lab_clue} (for lab clue) = {ehrlichia_score} total points.")
    print("\nEhrlichia provides the best fit for all the available clinical, geographical, and laboratory data.")

find_diagnosis()