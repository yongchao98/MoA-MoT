def calculate_max_handstand_score():
    """
    Calculates the maximum possible score for a MAG floor routine
    containing only handstand skills, based on the 2022 FIG Code of Points.
    """

    # 1. E-Score: Assumed to be perfect.
    e_score = 10.0

    # 2. D-Score Calculation
    # a. Skill Value: Top 10 highest-value, non-repeating handstand skills.
    # There are 6 D-value (0.4) and 4 C-value (0.3) skills selected.
    num_d_skills = 6
    d_skill_value = 0.4
    num_c_skills = 4
    c_skill_value = 0.3
    top_10_skills_total_value = (num_d_skills * d_skill_value) + (num_c_skills * c_skill_value)

    # b. Element Group (EG) Bonus: 0.5 for fulfilling EG I (non-acrobatic).
    eg_bonus = 0.5

    # c. Connection Value (CV): Maximized by performing three connections of 3 elements.
    # Each connection of 3 elements (B-value or higher) is worth 0.2.
    num_connections_of_3 = 3
    cv_per_connection_of_3 = 0.2
    connection_value = num_connections_of_3 * cv_per_connection_of_3

    # Total D-Score is the sum of the above components.
    d_score = top_10_skills_total_value + eg_bonus + connection_value

    # 3. Chief Judge (CJ) Deduction: A 0.5 deduction for no acrobatic elements.
    cj_deduction = 0.5

    # 4. Final Score Calculation
    final_score = d_score + e_score - cj_deduction

    # Print the detailed breakdown of the final score calculation.
    print("Calculating the Maximum Score for a Handstand-Only Floor Routine")
    print("-----------------------------------------------------------------")
    print(f"Execution Score (E-Score): {e_score:.1f} (Perfect Execution)")
    print("\nDifficulty Score (D-Score) Calculation:")
    print(f"  - Top 10 Skills Value: ({num_d_skills} * {d_skill_value:.1f}) + ({num_c_skills} * {c_skill_value:.1f}) = {top_10_skills_total_value:.1f}")
    print(f"  - Element Group Bonus (EG I): {eg_bonus:.1f}")
    print(f"  - Connection Value (3x connections of 3): {num_connections_of_3} * {cv_per_connection_of_3:.1f} = {connection_value:.1f}")
    print(f"  - Total D-Score: {top_10_skills_total_value:.1f} + {eg_bonus:.1f} + {connection_value:.1f} = {d_score:.1f}")
    print("\nChief Judge Deduction:")
    print(f"  - No Acrobatic Elements: {cj_deduction:.1f}")
    print("\n--- Final Score Calculation ---")
    print(f"Final Score = (D-Score) + (E-Score) - (CJ Deduction)")
    print(f"Final Score = {d_score:.1f} + {e_score:.1f} - {cj_deduction:.1f} = {final_score:.1f}")

calculate_max_handstand_score()
<<<14.2>>>