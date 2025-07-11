def calculate_max_handstand_score():
    """
    Calculates the highest possible score for a MAG floor routine
    consisting of only handstand skills, assuming perfect execution.
    """

    # 1. E-Score (Execution)
    # The prompt assumes perfect execution, which starts at 10.0.
    e_score = 10.0

    # 2. D-Score (Difficulty)
    # The D-Score has three components: Skill Value, Element Group Bonus, and Connection Value.

    # 2a. Skill Value (Top 10 skills)
    # We will construct a routine with the 10 highest-value, distinct handstand skills.
    # Based on the FIG Code of Points:
    # - 1x 'D' skill: Handstand with 1080° turn (D=0.4)
    # - 5x 'C' skills: e.g., 720° turn, various straight-arm presses from L-sit, felge to handstand (C=0.3)
    # - 4x 'B' skills: e.g., 360° turn, various bent-arm presses (B=0.2)
    
    d_skills_value = 1 * 0.4
    c_skills_value = 5 * 0.3
    b_skills_value = 4 * 0.2
    
    total_skill_value = d_skills_value + c_skills_value + b_skills_value

    # 2b. Element Group (EG) Bonus
    # The gymnast must fulfill 4 groups (I: Non-acro, II: Acro Fwd, III: Acro Bwd, IV: Acro Side/Other).
    # Since only Group I can be fulfilled, the 0.5 bonus is not awarded.
    eg_bonus = 0.0

    # 2c. Connection Value (CV)
    # We can connect our highest skills for a bonus.
    # Connection 1: Press to HS (C) + 1080° Turn (D) -> 0.2 CV
    # Connection 2: Press to HS (C) + 720° Turn (C) -> 0.2 CV
    # Connection 3: Felge to HS (C) + Press to HS (C) -> 0.2 CV
    connection_1 = 0.2
    connection_2 = 0.2
    connection_3 = 0.2
    total_connection_value = connection_1 + connection_2 + connection_3
    
    # Calculate the total D-Score
    d_score = total_skill_value + eg_bonus + total_connection_value
    
    # Calculate the Final Score
    final_score = d_score + e_score

    # Print the detailed breakdown
    print("--- Gymnastics Score Calculation ---")
    print("\nPart 1: Execution Score (E-Score)")
    print(f"Perfect Execution Score: {e_score}")
    
    print("\nPart 2: Difficulty Score (D-Score)")
    print("  Component A: Top 10 Skill Values")
    print(f"    1 'D' Skill  (1 * 0.4) = {d_skills_value:.1f}")
    print(f"    5 'C' Skills (5 * 0.3) = {c_skills_value:.1f}")
    print(f"    4 'B' Skills (4 * 0.2) = {b_skills_value:.1f}")
    print(f"    Total Skill Value = {d_skills_value:.1f} + {c_skills_value:.1f} + {b_skills_value:.1f} = {total_skill_value:.1f}")
    
    print("\n  Component B: Element Group (EG) Bonus")
    print(f"    EG Bonus (cannot be fulfilled): {eg_bonus}")

    print("\n  Component C: Connection Value (CV)")
    print(f"    Three C/D connections (0.2 + 0.2 + 0.2) = {total_connection_value:.1f}")

    print("\n  Total D-Score Calculation:")
    print(f"    D-Score = Skill Value + EG Bonus + CV")
    print(f"    D-Score = {total_skill_value:.1f} + {eg_bonus:.1f} + {total_connection_value:.1f} = {d_score:.1f}")
    
    print("\n--- Final Score ---")
    print(f"Final Score = D-Score + E-Score")
    print(f"Final Score = {d_score:.1f} + {e_score:.1f} = {final_score:.1f}")


if __name__ == '__main__':
    calculate_max_handstand_score()
    # The final answer is the numeric value of the calculated score.
    final_score = (1 * 0.4 + 5 * 0.3 + 4 * 0.2) + 0.0 + (0.2 + 0.2 + 0.2) + 10.0
    print(f"\n<<< {final_score:.1f} >>>")
