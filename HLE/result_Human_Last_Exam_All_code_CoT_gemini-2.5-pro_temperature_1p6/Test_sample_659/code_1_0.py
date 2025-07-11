def calculate_max_handstand_score():
    """
    Calculates the highest possible score for a MAG floor routine with only handstand skills,
    assuming perfect execution under the 2022-2024 FIG Code of Points.
    """
    
    # Base E-Score and Penalties
    e_score = 10.0
    no_dismount_penalty = 0.5

    print("--- Score Calculation Breakdown ---")
    print(f"\n1. E-Score (Execution):")
    print(f"Assuming perfect execution, the E-Score is {e_score:.1f}")

    # Step 1: Difficulty Value (DV) of the top 10 skills
    # Based on the FIG Code, we select the 10 highest-value unique handstand elements:
    # 2 'E' skills (0.5), 4 'D' skills (0.4), and 4 'C' skills (0.3).
    e_skills_count = 2
    e_value = 0.5
    d_skills_count = 4
    d_value = 0.4
    c_skills_count = 4
    c_value = 0.3

    dv_from_e = e_skills_count * e_value
    dv_from_d = d_skills_count * d_value
    dv_from_c = c_skills_count * c_value
    dv_total = dv_from_e + dv_from_d + dv_from_c
    
    print("\n2. Difficulty Value (DV) from Top 10 Skills:")
    print(f"   - {e_skills_count} 'E' skills x {e_value} = {dv_from_e:.1f}")
    print(f"   - {d_skills_count} 'D' skills x {d_value} = {dv_from_d:.1f}")
    print(f"   - {c_skills_count} 'C' skills x {c_value} = {dv_from_c:.1f}")
    print(f"   Total DV = {dv_from_e:.1f} + {dv_from_d:.1f} + {dv_from_c:.1f} = {dv_total:.1f}")

    # Step 2: Element Group (EG) Requirements Bonus
    # Only Group I (non-acrobatic) is fulfilled. Groups II, III (acrobatic) and IV (dismount) are not.
    eg_bonus = 0.5
    print("\n3. Element Group (EG) Bonus:")
    print(f"   Only Group I (Non-acrobatic) is fulfilled, awarding {eg_bonus:.1f} points.")
    
    # Step 3: Connection Value (CV) Bonus
    # Create a sequence of 10 skills to maximize connection values: E-E-D-D-D-D-C-C-C-C
    # This results in 9 connections.
    cv_ee = 0.2  # E -> E
    cv_ed = 0.2  # E -> D
    cv_dd = 0.2  # D -> D
    cv_dc = 0.1  # D -> C
    cv_cc = 0.1  # C -> C
    
    # The connections in the optimal sequence are:
    # E->E, E->D, D->D, D->D, D->D, D->C, C->C, C->C, C->C
    connections_calculation = [cv_ee, cv_ed, cv_dd, cv_dd, cv_dd, cv_dc, cv_cc, cv_cc, cv_cc]
    cv_total = sum(connections_calculation)
    
    print("\n4. Connection Value (CV) Bonus:")
    print(f"   Connecting the 10 skills in an optimal sequence yields 9 connections for a bonus:")
    print(f"   CV = {connections_calculation[0]} (E+E) + {connections_calculation[1]} (E+D) + {connections_calculation[2]} (D+D) + {connections_calculation[3]} (D+D) + {connections_calculation[4]} (D+D) + {connections_calculation[5]} (D+C) + {connections_calculation[6]} (C+C) + {connections_calculation[7]} (C+C) + {connections_calculation[8]} (C+C) = {cv_total:.1f}")

    # Step 4: Calculate final D-Score
    d_score = dv_total + eg_bonus + cv_total
    print("\n5. Total D-Score (Difficulty):")
    print(f"   D-Score = DV + EG + CV = {dv_total:.1f} + {eg_bonus:.1f} + {cv_total:.1f} = {d_score:.1f}")
    
    # Step 5: Final Score
    print("\n6. Final Score Calculation:")
    print(f"   A penalty is applied for lacking a valid dismount: -{no_dismount_penalty:.1f}")
    
    final_score = d_score + e_score - no_dismount_penalty

    print("\n--- Final Score ---")
    print(f"The final score is calculated as (D-Score + E-Score - Penalty).")
    print(f"Final Equation: {d_score:.1f} (D-Score) + {e_score:.1f} (E-Score) - {no_dismount_penalty:.1f} (Penalty) = {final_score:.1f}")

    return final_score

# Execute the calculation and print the final answer in the required format
final_score = calculate_max_handstand_score()
print(f"\n<<< {final_score} >>>")
