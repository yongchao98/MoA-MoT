def calculate_gymnastics_score():
    """
    Calculates the maximum possible score for a MAG floor routine
    using only handstand skills under the 2022 FIG Code of Points.
    """

    # 1. E-Score (Execution)
    # Assumes perfect execution as per the prompt.
    e_score = 10.0

    # 2. D-Score (Difficulty)
    # It is composed of Difficulty Value (DV), Composition Requirements (CR), and Connection Value (CV).

    # 2a. Difficulty Value (DV) - Top 10 handstand skills
    skill_values = {
        'E': 0.5, 'D': 0.4, 'C': 0.3, 'B': 0.2
    }
    # Top 10 valid skills: 1xE, 2xD, 4xC, 3xB
    dv_e = 1 * skill_values['E']
    dv_d = 2 * skill_values['D']
    dv_c = 4 * skill_values['C']
    dv_b = 3 * skill_values['B']
    difficulty_value = dv_e + dv_d + dv_c + dv_b

    # 2b. Composition Requirements (CR) - 2.5 total possible points
    # Only the requirement for a non-acrobatic element group is met.
    compositional_requirements = 0.5

    # 2c. Connection Value (CV)
    # Bonus for connecting difficult elements. D+C=0.1, D+D(or harder)=0.2
    # Optimal chain: C-D-E-D-C
    cv_cd1 = 0.1
    cv_de = 0.2
    cv_ed = 0.2
    cv_dc2 = 0.1
    connection_value = cv_cd1 + cv_de + cv_ed + cv_dc2

    # Calculate Total D-Score
    d_score = difficulty_value + compositional_requirements + connection_value

    # Calculate Final Score
    final_score = d_score + e_score

    # Print the breakdown of the calculation
    print("Score Calculation Breakdown:")
    print("==========================")
    print(f"E-Score (Execution): {e_score:.1f} (Perfect Execution)")
    print("\n--- D-Score Calculation ---")
    print("Difficulty Value (DV):")
    print(f"  - 1 'E' skill: 1 * {skill_values['E']} = {dv_e:.1f}")
    print(f"  - 2 'D' skills: 2 * {skill_values['D']} = {dv_d:.1f}")
    print(f"  - 4 'C' skills: 4 * {skill_values['C']} = {dv_c:.1f}")
    print(f"  - 3 'B' skills: 3 * {skill_values['B']} = {dv_b:.1f}")
    print(f"  Total DV = {dv_e:.1f} + {dv_d:.1f} + {dv_c:.1f} + {dv_b:.1f} = {difficulty_value:.1f}\n")

    print("Compositional Requirements (CR):")
    print(f"  - Only 1 of 5 requirements met = {compositional_requirements:.1f}\n")

    print("Connection Value (CV):")
    print(f"  - C+D connection: {cv_cd1:.1f}")
    print(f"  - D+E connection: {cv_de:.1f}")
    print(f"  - E+D connection: {cv_ed:.1f}")
    print(f"  - D+C connection: {cv_dc2:.1f}")
    print(f"  Total CV = {cv_cd1:.1f} + {cv_de:.1f} + {cv_ed:.1f} + {cv_dc2:.1f} = {connection_value:.1f}\n")

    print(f"Total D-Score = DV + CR + CV")
    print(f"Total D-Score = {difficulty_value:.1f} + {compositional_requirements:.1f} + {connection_value:.1f} = {d_score:.1f}\n")

    print("--- Final Score ---")
    print(f"Final Score = D-Score + E-Score")
    print(f"Final Score = {d_score:.1f} + {e_score:.1f} = {final_score:.1f}")

calculate_gymnastics_score()