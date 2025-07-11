import sys

# Define a function to encapsulate the diagnostic logic.
def run_diagnostic_analysis():
    """
    Analyzes clinical findings from a case study to determine the most likely diagnosis.
    """
    # Key findings from the case vignette
    location_of_pain = "left upper quadrant"
    referred_pain = "left shoulder" # This is Kehr's sign
    initial_hemoglobin = 11.7
    later_hemoglobin = 6.5
    was_polypectomy_performed = False
    hemodynamic_status = "unstable"

    # Define scores for each diagnosis, starting at 0.
    # We will add points based on how well the evidence fits.
    diagnostic_scores = {
        "Colonic perforation": 0,
        "Lower GI bleeding": 0,
        "Splenic laceration": 0,
        "Postpolypectomy syndrome": 0,
    }

    print("Step 1: Analyzing patient's symptoms and location of pain.")
    # Pain in the LUQ with referred pain to the left shoulder (Kehr's sign) is a classic sign of splenic injury.
    pain_location_score = 5
    if location_of_pain == "left upper quadrant" and referred_pain == "left shoulder":
        print(f"Finding: The patient has {location_of_pain} pain and {referred_pain} discomfort.")
        print(f"Conclusion: This combination is highly suggestive of splenic injury. Adding {pain_location_score} points to 'Splenic laceration'.\n")
        diagnostic_scores["Splenic laceration"] += pain_location_score

    print("Step 2: Analyzing hemodynamic status and blood loss.")
    # The drastic drop in hemoglobin indicates severe internal bleeding.
    hemoglobin_drop = initial_hemoglobin - later_hemoglobin
    hemorrhage_score = 4
    if hemoglobin_drop > 4 and hemodynamic_status == "unstable":
        print(f"Finding: Hemoglobin dropped significantly from {initial_hemoglobin} g/dL to {later_hemoglobin} g/dL, and the patient became hemodynamically unstable.")
        print(f"Conclusion: This indicates massive internal hemorrhage. Adding {hemorrhage_score} points to diagnoses involving major bleeding.\n")
        # Splenic laceration and colonic perforation can both cause massive bleeding.
        diagnostic_scores["Splenic laceration"] += hemorrhage_score
        diagnostic_scores["Colonic perforation"] += hemorrhage_score # A severe perforation could also bleed.
        diagnostic_scores["Lower GI bleeding"] += hemorrhage_score
    
    print("Step 3: Analyzing procedure history.")
    # The absence of polypectomy rules out postpolypectomy syndrome.
    # We will assign a very low score to effectively remove it from consideration.
    rule_out_score = -100
    if not was_polypectomy_performed:
        print("Finding: No polypectomy was performed during the colonoscopy.")
        print(f"Conclusion: This rules out 'Postpolypectomy syndrome'. Applying a large negative score of {rule_out_score}.\n")
        diagnostic_scores["Postpolypectomy syndrome"] += rule_out_score

    print("--- Final Diagnostic Calculation ---")
    print("The final scores for each diagnosis are calculated based on the evidence:")

    # Printing a representative "equation" for the most likely diagnosis
    # Final Score = (Pain Location Score) + (Hemorrhage Score)
    final_score_splenic_laceration = pain_location_score + hemorrhage_score
    print(f"Equation for Splenic Laceration: Score({pain_location_score}) + Score({hemorrhage_score}) = {final_score_splenic_laceration}\n")
    
    # Determine the most likely diagnosis by finding the max score.
    most_likely_diagnosis = max(diagnostic_scores, key=diagnostic_scores.get)

    print(f"Final scores: {diagnostic_scores}")
    print(f"The evidence most strongly supports a diagnosis of: {most_likely_diagnosis}")

# Entry point of the script
if __name__ == "__main__":
    original_stdout = sys.stdout  # Save a reference to the original standard output
    try:
        run_diagnostic_analysis()
        # The reasoning clearly points to 'Splenic laceration', which is option C.
        # According to the instructions, we should now output the final answer in the required format.
        print("\n<<<C>>>", file=original_stdout)

    except Exception as e:
        print(f"An error occurred: {e}", file=original_stdout)
