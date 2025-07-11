def diagnose_pelvic_mass():
    """
    Analyzes clinical findings to determine the most likely diagnosis
    for a 2-year-old boy with a pelvic mass.
    """

    # --- Patient's Clinical Findings ---
    # These are the key data points from the case.
    age = 2
    has_aniridia = True
    has_hypertension = True
    has_pelvic_mass = True
    has_developmental_delay = True # Represented by delayed speech

    # --- Diagnostic Scoring (for demonstration) ---
    # Assigning weights to findings to illustrate their importance.
    # The association between Aniridia and WAGR/Wilms Tumor is extremely high.
    aniridia_score = 10
    mass_at_age_2_score = 5
    hypertension_score = 3
    dev_delay_score = 4

    print("Analyzing the clinical case based on key findings...")
    print(f"Patient Age: {age} years")
    print(f"Finding 1: Aniridia -> Points to WAGR syndrome.")
    print(f"Finding 2: Pelvic Mass -> Consistent with Wilms Tumor.")
    print(f"Finding 3: Hypertension -> Consistent with Wilms Tumor (renin secretion).")
    print(f"Finding 4: Developmental Delay -> Consistent with WAGR syndrome.")
    print("-" * 30)

    # --- Constructing the Diagnostic 'Equation' ---
    # The problem asks to output the numbers in a final equation.
    # We will create a likelihood score for Nephroblastoma (Wilms Tumor).
    total_score = aniridia_score + mass_at_age_2_score + hypertension_score + dev_delay_score

    print("Calculating the likelihood score for Nephroblastoma (Wilms Tumor):")
    # This print statement fulfills the requirement to show each number in the final equation.
    print(f"Score = (Aniridia: {aniridia_score}) + (Pelvic Mass at Age {age}: {mass_at_age_2_score}) + (Hypertension: {hypertension_score}) + (Developmental Delay: {dev_delay_score})")
    print(f"Total Score = {total_score}")
    print("-" * 30)
    print("Conclusion:")
    print("The constellation of findings, especially the pathognomonic sign of aniridia, is definitive for WAGR syndrome.")
    print("The 'W' in WAGR stands for Wilms Tumor, which is a Nephroblastoma.")
    print("\nThe most likely diagnosis is Nephroblastoma.")


if __name__ == '__main__':
    diagnose_pelvic_mass()
<<<D>>>