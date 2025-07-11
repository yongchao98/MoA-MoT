def calculate_mag_floor_handstand_score():
    """
    Calculates the highest possible score for a MAG floor routine
    containing only handstand skills, assuming perfect execution.
    """

    # E-Score (Execution Score)
    # Assumes perfect execution, so no deductions from the base 10.0
    e_score = 10.0

    # D-Score (Difficulty Score) is composed of DV, EGR, and CV

    # 1. Difficulty Value (DV)
    # A gymnast can only count a maximum of 4 non-acrobatic (Group IV) elements.
    # The top 4 handstand skills are:
    # E-value skill (e.g., One-arm press to one-arm handstand): 0.5
    # E-value skill (e.g., 1080 turn on one arm): 0.5
    # E-value skill (e.g., Press from handstand thru planche to handstand): 0.5
    # D-value skill (e.g., Press to handstand from V-sit): 0.4
    top_4_handstand_skills = [0.5, 0.5, 0.5, 0.4]
    dv_score = sum(top_4_handstand_skills)

    # 2. Element Group Requirements (EGR)
    # The gymnast only fulfills the non-acrobatic element requirement.
    # EGR I (Non-Acro): 0.5
    # EGR II (Acro Forward): 0.0
    # EGR III (Acro Backward): 0.0
    # EGR IV (Dismount): 0.0
    egr_score = 0.5

    # 3. Connection Value (CV)
    # The highest connection bonus for two Group IV elements is for
    # connecting an E-value skill to another E-value skill.
    cv_score = 0.2

    # Calculate the total D-Score and Final Score
    d_score = dv_score + egr_score + cv_score
    total_score = d_score + e_score

    # Print the detailed breakdown of the score calculation
    print("Calculating the Maximum MAG Floor Score with Handstands Only:")
    print("-" * 60)
    print(f"Execution Score (E-Score): {e_score:.1f} (Perfect Execution)")
    print("\nDifficulty Score (D-Score):")
    print(f"  - Difficulty Value (DV): {top_4_handstand_skills[0]} + {top_4_handstand_skills[1]} + {top_4_handstand_skills[2]} + {top_4_handstand_skills[3]} = {dv_score:.1f}")
    print(f"  - Element Group Requirements (EGR): {egr_score:.1f} (Only Group I fulfilled)")
    print(f"  - Connection Value (CV): {cv_score:.1f} (Max bonus for non-acro connection)")
    print(f"  - Total D-Score = {dv_score:.1f} + {egr_score:.1f} + {cv_score:.1f} = {d_score:.1f}")
    print("-" * 60)
    print("\nFinal Score Calculation:")
    print("D-Score (DV + EGR + CV) + E-Score = Total Score")
    # Final print of the equation with all numbers
    print(f"({dv_score:.1f} + {egr_score:.1f} + {cv_score:.1f}) + {e_score:.1f} = {total_score:.1f}")


calculate_mag_floor_handstand_score()
<<<12.6>>>