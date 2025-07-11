def analyze_clinical_case():
    """
    This script analyzes the provided clinical case to determine the most likely diagnosis
    and the physical exam maneuver that would confirm it.
    """
    
    print("Step 1: Analyzing patient's key symptoms and signs.")
    # We assign numeric values to key findings to create a simple diagnostic score.
    # This is a conceptual model to illustrate the reasoning process.
    symptom_score = 1  # For L4-S1 radicular pain (sciatica).
    positional_pain_score = 1 # For pain intensified in specific positions (supine).
    imaging_score = 1  # For unremarkable X-ray, making bony pathology less likely.
    
    print("Patient presents with sciatica-like pain and positional discomfort.")
    print("Normal X-ray reduces likelihood of fracture or severe degenerative joint disease.")
    print("The total score from initial findings is: {} + {} + {} = 3".format(
        symptom_score, positional_pain_score, imaging_score
    ))
    
    print("\nStep 2: Evaluating the described physical exam maneuver.")
    print("The test is performed with the patient in the left decubitus position to test the right leg.")
    print("This position is ideal for isolating and testing the hip rotator muscles.")
    print("The primary function of the piriformis muscle is external rotation of the hip.")
    print("Piriformis Syndrome occurs when this muscle irritates the sciatic nerve.")
    
    print("\nStep 3: Determining the confirmatory action.")
    print("To confirm Piriformis Syndrome, a maneuver must stress the piriformis muscle to reproduce the pain.")
    print("Having the patient perform external rotation against resistance causes the piriformis to contract forcefully.")
    print("If this action reproduces the sciatica, it strongly confirms the diagnosis.")
    
    # We assign a higher weight to this highly specific test.
    specific_maneuver_score = 2
    
    print("\nStep 4: Calculating the final diagnostic likelihood score.")
    initial_score = symptom_score + positional_pain_score + imaging_score
    final_score = initial_score + specific_maneuver_score
    
    print("The confirmatory maneuver of resisted external rotation is highly specific.")
    # The final printout shows the "equation" as requested.
    print(f"Final Piriformis Syndrome Likelihood Score = {symptom_score} (symptoms) + {positional_pain_score} (positional) + {imaging_score} (imaging) + {specific_maneuver_score} (specific test) = {final_score}")

    print("\nConclusion: The action that will confirm the diagnosis of Piriformis Syndrome is external rotation against resistance.")
    
# Run the analysis
analyze_clinical_case()
<<<D>>>