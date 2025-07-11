def calculate_max_handstand_score():
    """
    Calculates the highest possible score for a MAG floor routine consisting only of handstand skills,
    based on the 2022 FIG Code of Points.
    """

    # --- Step 1: Define the optimal routine under the constraints ---
    # The gymnast is limited by:
    # 1. Max 3 Group I (strength) elements.
    # 2. Inability to perform a salto dismount.
    # We select the highest value skills that meet these rules. There are only 9 unique applicable skills.

    routine_skills = [
        # Top 3 Group I (Strength) Elements
        {'value': 0.5, 'name': "Manna press to handstand (E)"},
        {'value': 0.4, 'name': "V-sit press to handstand (D)"},
        {'value': 0.4, 'name': "Straddle L-sit, stoop, press to handstand (D)"},
        # Top Group III (Acrobatic Backward/Sideways) Handstand Elements
        {'value': 0.4, 'name': "Russian L-sit to handstand (D)"},
        {'value': 0.3, 'name': "Jump 1/1 turn to handstand (C)"},
        {'value': 0.3, 'name': "Yamawaki to handstand (C)"},
        {'value': 0.2, 'name': "Backward roll to handstand (B)"},
        {'value': 0.2, 'name': "Jump 1/2 turn to handstand (B)"},
        # Top Group II (Acrobatic Forward) Handstand Element
        {'value': 0.1, 'name': "Forward roll to handstand (A)"},
    ]

    skill_values = [skill['value'] for skill in routine_skills]

    # --- Step 2: Calculate the D-Score Components ---

    # Difficulty Value (DV) from the sum of all skill values
    total_dv = sum(skill_values)

    # Element Group (EG) Bonus
    # Group I, II, and III are fulfilled. Group IV (Dismount) is not.
    eg_bonus_g1 = 0.5
    eg_bonus_g2 = 0.5
    eg_bonus_g3 = 0.5
    total_eg_bonus = eg_bonus_g1 + eg_bonus_g2 + eg_bonus_g3

    # Connection Bonus (CB)
    # A +0.2 bonus is awarded for connecting two D-value or higher Group I strength elements.
    connection_bonus = 0.2

    # Deductions
    # A penalty is applied for not having a valid salto dismount.
    dismount_penalty = -0.5

    # Final D-Score
    d_score = total_dv + total_eg_bonus + connection_bonus + dismount_penalty

    # --- Step 3: Calculate the Final Score ---

    # E-Score is perfect by assumption
    e_score = 10.0

    # Final score is D-Score + E-Score
    final_score = d_score + e_score

    # --- Step 4: Output the Detailed Calculation ---

    print("To achieve the highest score with only handstand skills, the routine is constructed as follows:")
    print(f"\n1. The gymnast performs the 9 highest-value, unique handstand skills allowed under the rules.")
    print(f"2. The gymnast earns a perfect Execution Score (E-Score) of {e_score:.1f}.")
    print(f"3. The Difficulty Score (D-Score) is calculated based on skill values, bonuses, and penalties.")

    print("\n--- Final Score Calculation ---")

    # The final equation with each number explicitly shown
    skill_values_str = ' + '.join(f"{v:.1f}" for v in skill_values)
    eg_bonuses_str = f"{eg_bonus_g1:.1f} + {eg_bonus_g2:.1f} + {eg_bonus_g3:.1f}"

    print(f"\nFinal Score = (Difficulty Values) + (Element Group Bonuses) + (Connection Bonus) + (Penalty) + (E-Score)")
    print(f"Final Score = ({skill_values_str}) + ({eg_bonuses_str}) + {connection_bonus:.1f} {dismount_penalty:.1f} + {e_score:.1f}")
    print(f"Final Score = ({total_dv:.1f}) + ({total_eg_bonus:.1f}) + {connection_bonus:.1f} {dismount_penalty:.1f} + {e_score:.1f}")
    print(f"Final Score = {d_score:.1f} (D-Score) + {e_score:.1f} (E-Score)")
    print(f"Highest Possible Score = {final_score:.1f}")


calculate_max_handstand_score()