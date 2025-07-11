def calculate_max_handstand_score():
    """
    Calculates the maximum possible score for a MAG floor routine consisting
    only of handstand skills, assuming perfect execution under the 2022 FIG Code of Points.
    """

    # 1. E-Score (Execution)
    # Starts at 10.0. "Perfect execution" means no deductions.
    e_score = 10.0

    # 2. D-Score (Difficulty) = DV + EGR + CV

    # 2a. Difficulty Value (DV)
    # The top 10 most difficult, unique handstand skills from Group I.
    # Format: {"Skill Name": value}
    top_10_skills = {
        "V-sit to press to handstand with straight arms/body": 0.6,  # F-value
        "L-sit/straddle L-sit to press to HS w/ straight arms/body": 0.5, # E-value
        "Handstand, two full turns (720 deg) on one arm": 0.5,      # E-value
        "Press from stand to handstand w/ straight arms/body": 0.4,    # D-value
        "Handstand, one and a half turn (540 deg) on one arm": 0.4, # D-value
        "V-sit to press to handstand": 0.3,                           # C-value
        "Handstand, one full turn (360 deg) on one arm": 0.3,      # C-value
        "Straddle L-sit to press to handstand": 0.2,                  # B-value
        "Handstand with side split (straddle) legs, hold 2s": 0.2, # B-value
        "Handstand, hold 2s": 0.1                                     # A-value
    }
    
    skill_values = list(top_10_skills.values())
    dv_score = sum(skill_values)

    # 2b. Element Group Requirements (EGR)
    # A gymnast gets 0.5 points for each of 4 element groups fulfilled.
    # Handstands only fulfill Group I (Non-acrobatic elements).
    egr_score = 0.5

    # 2c. Connection Value (CV)
    # CV is only awarded for connecting acrobatic elements. Since the gymnast
    # only performs non-acrobatic skills, CV is 0.
    cv_score = 0.0

    # Calculate total D-Score
    d_score = dv_score + egr_score + cv_score

    # Calculate Final Score
    total_score = d_score + e_score

    # Print the explanation and breakdown
    print("Breakdown of the Maximum Floor Score with Handstand Skills Only:\n")
    print(f"E-Score (Execution): {e_score}")
    print("  - Based on the assumption of perfect execution.\n")
    
    print(f"D-Score (Difficulty): {d_score:.1f}")
    print("  - Calculated as DV + EGR + CV\n")

    print(f"  Difficulty Value (DV): {dv_score:.1f}")
    print("    - Sum of the 10 highest-valued unique handstand skills:")
    skill_values_str = " + ".join(map(str, sorted(skill_values, reverse=True)))
    print(f"    - DV Equation: {skill_values_str} = {dv_score:.1f}\n")

    print(f"  Element Group Requirements (EGR): {egr_score:.1f}")
    print("    - Only Group I (Non-Acrobatic) is fulfilled.\n")

    print(f"  Connection Value (CV): {cv_score:.1f}")
    print("    - Connection bonus is not awarded for non-acrobatic elements.\n")

    print("--- Final Score Calculation ---")
    # Final equation showing all numbers
    print(f"Total Score = E-Score + D-Score")
    print(f"{total_score:.1f} = {e_score:.1f} + {d_score:.1f} "
          f"(which is {dv_score:.1f} [DV] + {egr_score:.1f} [EGR] + {cv_score:.1f} [CV])")


calculate_max_handstand_score()
<<<14.0>>>