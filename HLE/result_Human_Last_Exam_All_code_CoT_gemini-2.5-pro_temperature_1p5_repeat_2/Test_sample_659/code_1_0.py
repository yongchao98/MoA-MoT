def calculate_max_handstand_score():
    """
    Calculates the maximum score for a MAG floor routine with only handstand skills,
    assuming perfect execution under the 2022 FIG Code of Points.
    """
    # 1. Execution Score (E-Score)
    e_score = 10.0

    # 2. Difficulty Score (D-Score) components

    # 2a. Difficulty Value (DV) from the top 9 skills
    # The highest value handstand skills are E=0.5, D=0.4, C=0.3.
    # To maximize DV, we take the 9 best available unique skills.
    # Assume 3 E-skills, 4 D-skills, and 2 C-skills.
    num_e_skills = 3
    val_e_skill = 0.5
    num_d_skills = 4
    val_d_skill = 0.4
    num_c_skills = 2
    val_c_skill = 0.3
    total_dv = (num_e_skills * val_e_skill) + (num_d_skills * val_d_skill) + (num_c_skills * val_c_skill)

    # 2b. Element Group Requirements (EGR)
    # Only EG I (Non-Acrobatic) is fulfilled.
    egr_I_value = 0.5
    total_egr = egr_I_value

    # 2c. Connection Value (CV)
    # Max of 4 connections are allowed. Connecting D+D or higher elements gives +0.2 each.
    num_connections = 4
    connection_bonus = 0.2
    total_cv = num_connections * connection_bonus

    # Calculate Total D-Score
    d_score = total_dv + total_egr + total_cv

    # Calculate Final Score
    final_score = d_score + e_score

    # Print the step-by-step calculation
    print("Calculating the highest possible score under the given conditions:\n")
    print(f"Final Score = Difficulty Score (D-Score) + Execution Score (E-Score)\n")
    print(f"Execution Score = {e_score:.1f} (due to 'perfect execution')\n")
    print(f"Difficulty Score = Difficulty Value (DV) + Element Group Requirements (EGR) + Connection Value (CV)\n")
    
    print(f"1. DV (top 9 skills): ({num_e_skills} * {val_e_skill}) + ({num_d_skills} * {val_d_skill}) + ({num_c_skills} * {val_c_skill}) = {total_dv:.1f}")
    print(f"2. EGR (only Group I fulfilled): {total_egr:.1f}")
    print(f"3. CV (4 connections): {num_connections} * {connection_bonus} = {total_cv:.1f}\n")

    print("--- Calculation ---")
    print(f"D-Score = {total_dv:.1f} + {total_egr:.1f} + {total_cv:.1f} = {d_score:.1f}")
    print(f"Final Score = {d_score:.1f} + {e_score:.1f} = {final_score:.1f}")

calculate_max_handstand_score()
<<<15.0>>>