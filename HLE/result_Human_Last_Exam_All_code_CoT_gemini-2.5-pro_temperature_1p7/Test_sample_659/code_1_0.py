def calculate_max_handstand_score():
    """
    Calculates the maximum score for a MAG floor routine with only handstand skills,
    assuming perfect execution under the 2022 FIG Code of Points.
    """

    # 1. E-Score (Execution)
    # The prompt assumes perfect execution, which starts at 10.0.
    e_score = 10.0

    # 2. D-Score (Difficulty) Components

    # 2a. Skill Difficulty (Value of the top 10 elements)
    # We select the 10 highest-value handstand skills available in the Code of Points.
    # The highest values are F(0.6), E(0.5), and D(0.4).
    # Optimal selection: 2x F-skills, 3x E-skills, 5x D-skills.
    skill_values = [0.6, 0.6, 0.5, 0.5, 0.5, 0.4, 0.4, 0.4, 0.4, 0.4]
    skill_difficulty_total = sum(skill_values)

    # 2b. Element Group (EG) Bonus
    # A gymnast gets 0.5 for each fulfilled requirement (Groups I, II, III).
    # All handstand skills are non-acrobatic (Group I). Groups II and III are not fulfilled.
    eg_bonus = 0.5

    # 2c. Connection Value (CV)
    # The CoP for floor exercise does not award CV for connecting high-value non-acrobatic skills.
    cv_bonus = 0.0

    # Calculate total D-Score
    d_score = skill_difficulty_total + eg_bonus + cv_bonus

    # 3. Final Score
    final_score = d_score + e_score

    # --- Output Results ---
    print("Calculating the Maximum Score for a Handstand-Only MAG Floor Routine")
    print("======================================================================")

    print(f"\nExecution Score (E-Score) = {e_score:.1f}")

    print("\n--- Difficulty Score (D-Score) Breakdown ---")
    
    # Print the equation for skill difficulty
    skill_values_str = " + ".join(map(str, skill_values))
    print(f"1. Skill Values (Top 10): {skill_values_str} = {skill_difficulty_total:.1f}")

    # Print the equation for EG bonus
    print(f"2. Element Group Bonus: {eg_bonus:.1f} (only Group I fulfilled)")

    # Print the equation for CV
    print(f"3. Connection Value (CV): {cv_bonus:.1f}")

    # Print the D-Score calculation
    print("--------------------------------------------------")
    d_score_equation = f"{skill_difficulty_total:.1f} (Skills) + {eg_bonus:.1f} (EG) + {cv_bonus:.1f} (CV)"
    print(f"Total D-Score = {d_score_equation} = {d_score:.1f}")
    print("--------------------------------------------------")

    print("\n--- Final Score Calculation ---")
    final_score_equation = f"{d_score:.1f} (D-Score) + {e_score:.1f} (E-Score)"
    print(f"Final Score = {final_score_equation} = {final_score:.1f}")
    print("======================================================================")


if __name__ == '__main__':
    calculate_max_handstand_score()