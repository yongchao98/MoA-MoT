def calculate_max_handstand_score():
    """
    Calculates the highest possible score for a MAG floor routine consisting
    only of handstand skills, based on the 2022 FIG Code of Points.
    """

    # 1. Top 8 Skill Values
    # A senior routine's difficulty is based on the top 8 skills.
    # The highest value handstand skills on floor are D(0.4), C(0.3), and B(0.2).
    # The gymnast can perform at least 3 unique D-skills, 3 unique C-skills,
    # and 2 unique B-skills.
    top_8_skills = [0.4, 0.4, 0.4, 0.3, 0.3, 0.3, 0.2, 0.2]
    skill_value_sum = sum(top_8_skills)

    # 2. Connection Bonus (CV)
    # The highest connection bonus for Group I elements is by connecting three
    # D-rated skills (D + D + D).
    connection_bonus = 0.3

    # 3. Penalties (Deductions from D-Score)
    # The routine must satisfy Element Group Requirements (EGR).
    # Since only handstands (Group I) are performed, there are penalties.
    penalty_group_II = 0.5  # For missing acrobatic elements forward
    penalty_group_III = 0.5 # For missing acrobatic elements backward
    penalty_no_dismount = 0.5 # For absence of a required acrobatic dismount
    
    penalties_list = [penalty_group_II, penalty_group_III, penalty_no_dismount]
    penalty_sum = sum(penalties_list)

    # 4. D-Score Calculation
    # D-Score = (Sum of Top 8 Skills) + (Connection Bonus) - (Penalties)
    d_score = skill_value_sum + connection_bonus - penalty_sum

    # 5. E-Score
    # The problem assumes perfect execution.
    e_score = 10.0

    # 6. Final Score
    final_score = d_score + e_score
    
    # 7. Print the final breakdown as a single equation
    skills_str = " + ".join(map(str, top_8_skills))
    penalties_str = " + ".join(map(str, penalties_list))
    
    print("This script calculates the maximum possible score for a MAG floor routine with only handstand skills.")
    print("\nThe final score is the sum of the Difficulty Score (D-Score) and the Execution Score (E-Score).\n")
    print("---")
    print("D-Score Calculation:")
    print(f"  Value of Top 8 Skills: {skills_str} = {skill_value_sum:.1f}")
    print(f"  Connection Bonus (CV): {connection_bonus:.1f}")
    print(f"  Penalties (Groups + Dismount): -({penalties_str}) = -{penalty_sum:.1f}")
    print(f"Total D-Score = {skill_value_sum:.1f} + {connection_bonus:.1f} - {penalty_sum:.1f} = {d_score:.1f}")
    print("\nE-Score (assuming perfect execution): 10.0")
    print("---\n")
    
    # Final combined equation
    print("Final Score Equation:")
    print(f"(({skills_str}) + {connection_bonus:.1f} - ({penalties_str})) + {e_score:.1f} = {final_score:.1f}")

calculate_max_handstand_score()
