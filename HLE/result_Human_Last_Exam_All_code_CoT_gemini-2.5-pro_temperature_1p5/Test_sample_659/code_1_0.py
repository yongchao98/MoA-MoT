def calculate_max_handstand_score():
    """
    Calculates the maximum possible score for a MAG floor routine
    consisting of only handstand skills, assuming perfect execution and
    that connection values for non-acrobatic skills are available.
    """

    # --- 1. D-Score: Skill Values ---
    # We select the top 10 available handstand skills from the Code of Points.
    # 1 x D-value skill (One-arm handstand)
    # 2 x C-value skills (e.g., Manna press to handstand)
    # 4 x B-value skills
    # 3 x A-value skills (to reach the 10-skill limit)
    d_skills = 1 * 0.4
    c_skills = 2 * 0.3
    b_skills = 4 * 0.2
    a_skills = 3 * 0.1
    total_skill_value = d_skills + c_skills + b_skills + a_skills

    # --- 2. D-Score: Element Group (EG) Bonus ---
    # Handstand skills are non-acrobatic and fall into Element Group I.
    # No acrobatic skills are performed, so only one of the four groups is fulfilled.
    eg_bonus = 0.5

    # --- 3. D-Score: Connection Value (CV) Bonus ---
    # Based on the prompt's assumption, we apply acrobatic connection rules.
    # To maximize CV, we arrange the 10 skills (D,C,C,B,B,B,B,A,A,A) in an optimal sequence.
    # A high-yield sequence is B-D-C-C-B-B-B-A-A-A, creating 9 connections:
    # 1. B -> D (B+D): 0.2 CV
    # 2. D -> C (C+D): 0.2 CV
    # 3. C -> C (C+C): 0.1 CV
    # 4. C -> B (B+C): 0.1 CV
    # 5. B -> B (B+B): 0.1 CV
    # 6. B -> B (B+B): 0.1 CV
    # 7. B -> A (A+B): 0.1 CV
    # 8. A -> A (A+A): 0.0 CV
    # 9. A -> A (A+A): 0.0 CV
    connection_values = [0.2, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.0, 0.0]
    cv_bonus = sum(connection_values)

    # --- 4. Total D-Score ---
    d_score = total_skill_value + eg_bonus + cv_bonus

    # --- 5. E-Score ---
    # The problem assumes perfect execution.
    e_score = 10.0

    # --- Final Score ---
    final_score = d_score + e_score

    # --- Print the explanation and final calculation ---
    print("Calculating the Maximum MAG Floor Score with Only Handstands:")
    print("-" * 60)
    print("1. D-Score Calculation:")
    print(f"   - Skill Values (Top 10): {total_skill_value:.1f}")
    print(f"   - Element Group Bonus (EG I only): {eg_bonus:.1f}")
    print(f"   - Connection Value Bonus (special condition): {cv_bonus:.1f}")
    print(f"   Total D-Score = {total_skill_value:.1f} + {eg_bonus:.1f} + {cv_bonus:.1f} = {d_score:.1f}\n")

    print("2. E-Score Calculation:")
    print(f"   - Assumed Perfect Execution: {e_score:.1f}\n")

    print("3. Final Score Calculation:")
    print("   D-Score + E-Score = Final Score")
    # The final equation with each number outputted
    print(f"   {d_score:.1f} + {e_score:.1f} = {final_score:.1f}\n")

    print(f"The highest achievable score is {final_score:.1f}.")

if __name__ == '__main__':
    calculate_max_handstand_score()
