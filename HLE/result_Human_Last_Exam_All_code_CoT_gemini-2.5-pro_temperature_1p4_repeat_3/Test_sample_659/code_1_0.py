def calculate_max_handstand_score():
    """
    Calculates the highest possible score for a MAG floor routine
    containing only handstand skills, assuming perfect execution under the
    2022-2024 FIG Code of Points.
    """

    # 1. E-Score (Execution)
    # Starts at 10.0 for perfect execution.
    e_score = 10.0

    # 2. D-Score (Difficulty) Components

    # 2a. Element Group (EG) Requirements
    # A routine needs skills from Group I (non-acrobatic), II (fwd acrobatic),
    # and III (bwd acrobatic).
    # Only Group I is fulfilled with handstand skills. 0.5 points are awarded.
    eg_points = 0.5

    # 2b. Difficulty Value (DV) from the top 10 skills
    # We find the 10 highest-value unique handstand-related skills in the Code of Points.
    # Highest value skills are presses to handstand.
    num_d_skills = 3
    d_skill_value = 0.4
    num_c_skills = 7
    c_skill_value = 0.3

    dv_from_d_skills = num_d_skills * d_skill_value
    dv_from_c_skills = num_c_skills * c_skill_value
    total_dv = dv_from_d_skills + dv_from_c_skills

    # 2c. Connection Value (CV)
    # A bonus of 0.2 is awarded for connecting two D-value strength elements
    # performed from different positions (e.g., press from V-sit connected to a
    # press from straddle-sit).
    cv_points = 0.2

    # Calculate Total D-Score
    d_score = total_dv + eg_points + cv_points

    # Calculate Final Score
    final_score = e_score + d_score

    # --- Print the step-by-step calculation ---
    print("--- MAG Floor Score Calculation (Handstand-Only Routine) ---\n")
    print("The final score is the sum of the E-Score and the D-Score.\n")

    print(f"E-Score (Execution):")
    print(f"Based on perfect execution, the E-Score is {e_score}\n")

    print("D-Score (Difficulty):")
    print("The D-Score is the sum of Difficulty Value (DV), Element Groups (EG), and Connection Value (CV).\n")

    print("1. Difficulty Value (DV):")
    print(f"The top 10 skills include:")
    print(f"- {num_d_skills} D-value skills worth {d_skill_value} each = {dv_from_d_skills:.1f} points")
    print(f"- {num_c_skills} C-value skills worth {c_skill_value} each = {dv_from_c_skills:.1f} points")
    print(f"Total DV = {dv_from_d_skills:.1f} + {dv_from_c_skills:.1f} = {total_dv:.1f} points\n")

    print("2. Element Groups (EG):")
    print(f"Only Group I (non-acrobatic) is fulfilled = {eg_points} points\n")

    print("3. Connection Value (CV):")
    print(f"Connecting two D-value strength elements = {cv_points} points\n")

    print("Total D-Score Calculation:")
    print(f"D-Score = DV ({total_dv:.1f}) + EG ({eg_points}) + CV ({cv_points})")
    print(f"D-Score = {d_score:.1f}\n")

    print("--- Final Score Calculation ---")
    print(f"Final Score = E-Score + D-Score")
    print(f"Final Score = {e_score:.1f} + {d_score:.1f} = {final_score:.1f}")


if __name__ == "__main__":
    calculate_max_handstand_score()
