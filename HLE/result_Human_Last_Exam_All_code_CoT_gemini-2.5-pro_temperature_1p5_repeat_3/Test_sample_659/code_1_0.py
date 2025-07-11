import sys
# Redirect stdout to a file to suppress output in this block
original_stdout = sys.stdout
sys.stdout = open('output.txt', 'w')
# --- Your code that you don't want to output ---
print("This output is suppressed.")
# --- Restore stdout ---
sys.stdout.close()
sys.stdout = original_stdout

def solve_gymnastics_score():
    """
    Calculates the highest possible score for a MAG floor routine
    containing only handstand skills, based on the 2022 FIG Code of Points.
    """

    # --- E-Score (Execution) ---
    # Assuming perfect execution as per the prompt.
    e_score = 10.0

    # --- D-Score (Difficulty) Components ---

    # 1. Top 10 Skill Values
    # Based on the 2022 FIG MAG Code of Points, these are the 10 highest-valued
    # unique non-acrobatic strength/handstand elements for floor.
    # G=0.7, F=0.6, E=0.5, D=0.4
    skill_values = [0.7, 0.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.4, 0.4]
    d_score_from_skills = sum(skill_values)

    # 2. Element Group (EG) Bonus
    # Only Element Group I (Non-Acrobatic Skills) is fulfilled.
    # Groups II (Fwd Acro) and III (Bwd Acro) are not.
    eg_bonus = 0.5

    # 3. Connection Value (CV) Bonus
    # The CoP awards 0.2 points for connecting two D-value (0.4) or higher strength elements.
    # All 10 chosen skills are D-value or higher.
    # A sequence of 10 elements allows for 9 direct connections.
    num_elements = len(skill_values)
    num_connections = num_elements - 1
    cv_per_connection = 0.2
    cv_total_bonus = num_connections * cv_per_connection

    # --- Final Score Calculation ---
    d_score = d_score_from_skills + eg_bonus + cv_total_bonus
    final_score = d_score + e_score

    # --- Print the detailed breakdown ---
    print("To determine the highest score, we calculate the D-Score and add the E-Score.")
    print("-------------------------------------------------------------------------")

    # D-Score Breakdown
    print("1. D-Score (Difficulty) Calculation")
    print("\nD-Score = (Sum of Top 10 Skills) + (Element Group Bonus) + (Connection Bonus)\n")

    skill_values_str = " + ".join(map(str, sorted(skill_values, reverse=True)))
    print(f"  A. Sum of Skills: The top 10 handstand skills give a value of:")
    print(f"     {skill_values_str} = {d_score_from_skills:.1f}\n")

    print(f"  B. Element Group Bonus: The routine only fulfills Group I (Non-Acrobatic).")
    print(f"     Bonus = {eg_bonus:.1f}\n")

    print(f"  C. Connection Bonus: Connecting all {num_elements} D+ skills yields {num_connections} connections.")
    print(f"     {num_connections} connections * {cv_per_connection} points = {cv_total_bonus:.1f}\n")

    print(f"  Total D-Score = {d_score_from_skills:.1f} (Skills) + {eg_bonus:.1f} (EG) + {cv_total_bonus:.1f} (CV)")
    print(f"  Total D-Score = {d_score:.1f}\n")

    print("-------------------------------------------------------------------------")

    # E-Score and Final Score
    print("2. E-Score (Execution) Calculation")
    print(f"With perfect execution, the E-Score is {e_score:.1f}.\n")

    print("-------------------------------------------------------------------------")
    print("3. Final Score Calculation")
    print("Final Score = D-Score + E-Score")
    print(f"Final Score = {d_score:.1f} + {e_score:.1f}")
    print(f"\nHighest Possible Score = {final_score:.1f}")


solve_gymnastics_score()
final_score = 17.4
print(f'<<<{final_score}>>>')