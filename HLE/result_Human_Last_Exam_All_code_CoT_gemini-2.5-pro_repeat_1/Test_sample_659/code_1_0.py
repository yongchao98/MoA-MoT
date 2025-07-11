def calculate_max_handstand_score():
    """
    Calculates the highest possible score for a MAG floor routine
    containing only handstand skills, assuming perfect execution.
    """

    # Part 1: E-Score (Execution Score)
    # The prompt assumes perfect execution, so the E-Score is 10.0.
    e_score = 10.0

    # Part 2: D-Score (Difficulty Score)

    # Step 2a: Sum of Difficulty Values (DV)
    # Handstand skills are non-acrobatic (Group I). We select the 10
    # highest-value unique skills. The dismount is one of these 10 skills.
    # F=0.6, E=0.5, D=0.4, C=0.3, B=0.2, A=0.1
    top_10_skills = {
        "Handstand with 720° (double) turn on one arm": 0.6,  # F-value
        "Handstand with 540° (1.5) turn on one arm": 0.5,     # E-value
        "V-sit, press to handstand with straight arms": 0.4, # D-value
        "Manna hold for 2 seconds": 0.4,                     # D-value
        "Planche hold for 2 seconds": 0.4,                   # D-value
        "Handstand with 360° (full) turn on one arm": 0.4,   # D-value
        "L-sit, press to handstand with straight arms": 0.3, # C-value
        "V-sit hold for 2 seconds": 0.3,                     # C-value
        "Straddle planche hold for 2 seconds": 0.3,          # C-value
        "Handstand, fall to Japanese HS to stand (Dismount)": 0.2, # B-value
    }
    sum_dv = sum(top_10_skills.values())

    # Step 2b: Connection Value (CV)
    # Connection Value is not awarded for linking two non-acrobatic (Group I) skills.
    cv = 0.0

    # Step 2c: Composition Requirement (CR) Deductions
    # A 0.5 deduction is applied for each missing requirement.
    # CR1: Fulfilled (routine is all Group I). Deduction = 0.0
    # CR2: Not fulfilled (needs Group II forward acro). Deduction = 0.5
    # CR3: Not fulfilled (needs Group III backward acro). Deduction = 0.5
    # CR4: Not fulfilled (dismount must be 'C' or higher, but is 'B'). Deduction = 0.5
    cr_deductions = 0.5 + 0.5 + 0.5

    # Step 2d: Final D-Score Calculation
    d_score = sum_dv + cv - cr_deductions

    # Part 3: Final Score Calculation
    final_score = d_score + e_score

    # Print the detailed breakdown
    print("--- Calculating Maximum Score for a Handstand-Only Floor Routine ---")
    print("\nPart 1: E-Score (Execution)")
    print(f"Assuming perfect execution, the base E-Score is {e_score:.1f}")

    print("\nPart 2: D-Score (Difficulty)")
    print("\nStep 2a: Sum of Difficulty Values (DV) of the Top 10 Skills:")
    skill_dv_list = [v for v in top_10_skills.values()]
    print("The 10 highest skills are valued at: "
          f"{skill_dv_list[0]:.1f}, {skill_dv_list[1]:.1f}, {skill_dv_list[2]:.1f}, "
          f"{skill_dv_list[3]:.1f}, {skill_dv_list[4]:.1f}, {skill_dv_list[5]:.1f}, "
          f"{skill_dv_list[6]:.1f}, {skill_dv_list[7]:.1f}, {skill_dv_list[8]:.1f}, "
          f"and {skill_dv_list[9]:.1f}")
    
    dv_sum_str = " + ".join(map(str, skill_dv_list))
    print(f"Sum of DV = {dv_sum_str} = {sum_dv:.1f}")


    print("\nStep 2b: Connection Value (CV)")
    print("No connections are awarded for non-acrobatic skills.")
    print(f"CV = {cv:.1f}")

    print("\nStep 2c: Composition Requirement (CR) Deductions")
    print("The routine is missing 3 of the 4 requirements:")
    print("- Missing forward acrobatic element: -0.5")
    print("- Missing backward acrobatic element: -0.5")
    print("- Dismount value is too low ('B' instead of 'C'+): -0.5")
    print(f"Total CR Deductions = 0.5 + 0.5 + 0.5 = {cr_deductions:.1f}")

    print("\nStep 2d: Final D-Score")
    print(f"D-Score = (Sum of DV) + (CV) - (CR Deductions)")
    print(f"D-Score = {sum_dv:.1f} + {cv:.1f} - {cr_deductions:.1f} = {d_score:.1f}")

    print("\n--- Part 3: Final Score ---")
    print("Final Score = D-Score + E-Score")
    print(f"Final Score = {d_score:.1f} + {e_score:.1f} = {final_score:.1f}")


calculate_max_handstand_score()
<<<12.3>>>