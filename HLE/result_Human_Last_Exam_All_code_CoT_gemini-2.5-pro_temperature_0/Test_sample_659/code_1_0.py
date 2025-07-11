def calculate_max_handstand_score():
    """
    Calculates the maximum possible score for a MAG floor routine
    consisting of only handstand skills, based on the 2022 FIG Code of Points.
    """
    # All available handstand skills on Floor with their difficulty values
    # (Element Name, Value)
    all_skills = [
        ("Handstand hop with 1080° turn (F)", 0.6),
        ("Press from Manna to handstand (E)", 0.5),
        ("Press from planche to handstand (E)", 0.5),
        ("Handstand hop with 900° turn (E)", 0.5),
        ("Press from V-sit to handstand (D)", 0.4),
        ("Press from straddle planche to handstand (D)", 0.4),
        ("Handstand with 720° turn on one arm (D)", 0.4),
        ("Handstand hop with 720° turn (D)", 0.4),
        ("Press from straddle L-sit to handstand (C)", 0.3),
        ("Handstand with 360° turn on one arm (C)", 0.3),
        ("Handstand hop with 540° turn (C)", 0.3),
        ("Press to handstand with straight arms (B)", 0.2),
        ("Handstand hop with 360° turn (B)", 0.2),
        ("Press to handstand (A)", 0.1),
        ("Handstand, hold 2 sec. (A)", 0.1),
        ("Handstand hop with 180° turn (A)", 0.1)
    ]

    # Sort skills by value in descending order and select the top 10
    top_10_skills = sorted(all_skills, key=lambda x: x[1], reverse=True)[:10]

    # --- D-Score Calculation ---
    print("Calculating the Difficulty Score (D-Score):")
    print("-" * 40)

    # 1. Sum of top 10 skill values
    skill_value_sum = sum(skill[1] for skill in top_10_skills)
    print("1. Top 10 Skill Values:")
    skill_values_str = []
    for name, value in top_10_skills:
        print(f"   - {name}: {value:.1f}")
        skill_values_str.append(f"{value:.1f}")
    
    print(f"\n   Sum of Skill Values = {' + '.join(skill_values_str)} = {skill_value_sum:.1f}")
    print("-" * 40)

    # 2. Element Group Requirements (EGR)
    egr_i_bonus = 0.5  # Fulfilled (all skills are non-acrobatic)
    egr_dismount_bonus = 0.5 # Fulfilled (can use a D+ non-acro skill)
    egr_total_bonus = egr_i_bonus + egr_dismount_bonus
    print("2. Element Group Requirement (EGR) Bonus:")
    print(f"   - Group I (Non-Acrobatic): Fulfilled = {egr_i_bonus:.1f}")
    print(f"   - Group II (Acrobatic Fwd): Not Fulfilled = 0.0")
    print(f"   - Group III (Acrobatic Bwd): Not Fulfilled = 0.0")
    print(f"   - Dismount: Fulfilled = {egr_dismount_bonus:.1f}")
    print(f"\n   Total EGR Bonus = {egr_i_bonus:.1f} + {egr_dismount_bonus:.1f} = {egr_total_bonus:.1f}")
    print("-" * 40)

    # 3. Connection Value (CV)
    connection_value = 0.0
    print("3. Connection Value (CV):")
    print("   - No CV for Non-Acrobatic -> Non-Acrobatic connections on Floor.")
    print(f"   Total CV = {connection_value:.1f}")
    print("-" * 40)

    # Calculate D-Score
    d_score = skill_value_sum + egr_total_bonus + connection_value
    print(f"Total D-Score = (Skill Values) + (EGR) + (CV)")
    print(f"Total D-Score = {skill_value_sum:.1f} + {egr_total_bonus:.1f} + {connection_value:.1f} = {d_score:.1f}")
    print("=" * 40)

    # --- E-Score and Final Score ---
    e_score = 10.0
    print("Calculating the Final Score:")
    print("-" * 40)
    print(f"Execution Score (E-Score) for perfect execution = {e_score:.1f}")
    
    final_score = d_score + e_score
    print(f"\nFinal Score = D-Score + E-Score")
    print(f"Final Score = {d_score:.1f} + {e_score:.1f} = {final_score:.1f}")
    print("=" * 40)

    return final_score

if __name__ == '__main__':
    highest_score = calculate_max_handstand_score()
    # The final answer is printed in the required format below.
    # The function above provides the detailed breakdown.
    print(f"\nThe highest possible score is {highest_score:.1f}.")
    print(f"<<<{highest_score:.1f}>>>")
