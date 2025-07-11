def calculate_max_handstand_score():
    """
    Calculates the highest possible score for a MAG floor routine
    containing only handstand skills, assuming perfect execution.
    """

    # 1. E-Score (Execution)
    # Assumes perfect execution as per the problem description.
    e_score = 10.0

    # 2. D-Score (Difficulty) Calculation

    # 2a. Difficulty Value (DV) from skills
    # The MAG code limits a floor routine to a maximum of 3 Element Group I (non-acrobatic) elements.
    # We select the three highest-valued handstand skills.
    # Skill 1: V-cross to handstand (E-value)
    skill_1_val = 0.5
    # Skill 2: Manna to handstand (E-value)
    skill_2_val = 0.5
    # Skill 3: Felge to handstand (D-value)
    skill_3_val = 0.4
    
    dv = skill_1_val + skill_2_val + skill_3_val

    # 2b. Element Group (EG) Requirements
    # There are 4 EGs worth 0.5 each.
    # EG I (Non-acrobatic): Fulfilled by the handstand skills.
    # EG II (Acrobatic Forward): Not fulfilled.
    # EG III (Acrobatic Backward): Not fulfilled.
    # EG IV (Dismount): Not fulfilled.
    eg_bonus = 0.5

    # 2c. Connection Value (CV)
    # Connecting the two E-value skills (D+D or higher) yields a 0.2 bonus.
    cv_bonus = 0.2

    # Total D-Score
    d_score = dv + eg_bonus + cv_bonus

    # 3. Neutral Deductions
    # A 0.5 deduction is applied for not performing a dismount.
    dismount_deduction = 0.5

    # 4. Final Score
    final_score = d_score + e_score - dismount_deduction

    # Output the results
    print("This script calculates the maximum MAG Floor score with only handstand skills.")
    print("\n--- Score Breakdown ---")
    print(f"Difficulty Value (DV): {skill_1_val} + {skill_2_val} + {skill_3_val} = {dv:.1f}")
    print(f"Element Group Bonus (EG): {eg_bonus:.1f} (Only Group I is fulfilled)")
    print(f"Connection Value Bonus (CV): {cv_bonus:.1f} (Connecting two E-value skills)")
    print(f"Total D-Score: {dv:.1f} (DV) + {eg_bonus:.1f} (EG) + {cv_bonus:.1f} (CV) = {d_score:.1f}")
    print(f"\nExecution Score (E-Score): {e_score:.1f} (Perfect execution)")
    print(f"Neutral Deductions: {dismount_deduction:.1f} (No dismount performed)")

    print("\n--- Final Score Calculation ---")
    # The final print statement showing the full equation as requested.
    print(f"{d_score:.1f} + {e_score:.1f} - {dismount_deduction:.1f} = {final_score:.1f}")


if __name__ == "__main__":
    calculate_max_handstand_score()
    # The final answer is derived from the calculation
    # D-Score(2.1) + E-Score(10.0) - Deduction(0.5) = 11.6
    print("\n<<<11.6>>>")