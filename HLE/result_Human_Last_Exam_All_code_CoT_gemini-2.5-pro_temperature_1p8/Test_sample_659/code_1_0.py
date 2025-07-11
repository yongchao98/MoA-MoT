def calculate_max_handstand_score():
    """
    Calculates the highest possible score for a MAG floor routine
    consisting of only handstand skills, assuming perfect execution.
    """

    # 1. Execution Score (E-Score)
    # Assumed to be perfect.
    e_score = 10.0

    # 2. Difficulty Score (D-Score) Components

    # 2a. Skill Value (top 10 skills)
    # We assume the gymnast performs the 10 highest-value distinct handstand skills.
    # A plausible distribution from the Code of Points:
    # 2 'D' skills (e.g., handstand with 720 turn), 4 'C' skills (e.g., various presses), 4 'B' skills.
    num_d_skills = 2
    val_d_skill = 0.4
    num_c_skills = 4
    val_c_skill = 0.3
    num_b_skills = 4
    val_b_skill = 0.2

    skill_value = (num_d_skills * val_d_skill) + \
                  (num_c_skills * val_c_skill) + \
                  (num_b_skills * val_b_skill)

    # 2b. Element Group Requirements (EGR)
    # The gymnast only performs Group I elements (non-acrobatic/handstands).
    # He misses Groups II, III, and IV. So, he only gets points for Group I.
    egr_points = 0.5

    # 2c. Connection Value (CV)
    # Assuming a connection between two high-value strength elements (e.g., C + D).
    connection_bonus = 0.2

    # Calculate total D-Score
    d_score = skill_value + egr_points + connection_bonus

    # Calculate Final Score
    final_score = d_score + e_score

    # Print the breakdown of the calculation
    print("Calculating the Maximum Score:")
    print("-" * 30)
    print(f"Execution Score (E-Score): {e_score}")
    print("-" * 30)
    print("Difficulty Score (D-Score) Calculation:")
    print(f"  Skill Value: ({num_d_skills} * {val_d_skill}) + ({num_c_skills} * {val_c_skill}) + ({num_b_skills} * {val_b_skill}) = {skill_value:.1f}")
    print(f"  Element Group Points: {egr_points} (Only Group I fulfilled)")
    print(f"  Connection Bonus: {connection_bonus}")
    print(f"Total D-Score Equation: {skill_value:.1f} + {egr_points} + {connection_bonus} = {d_score:.1f}")
    print("-" * 30)
    print(f"Final Score Equation (D-Score + E-Score):")
    print(f"{d_score:.1f} + {e_score} = {final_score:.1f}")
    print("-" * 30)

calculate_max_handstand_score()
<<<13.5>>>