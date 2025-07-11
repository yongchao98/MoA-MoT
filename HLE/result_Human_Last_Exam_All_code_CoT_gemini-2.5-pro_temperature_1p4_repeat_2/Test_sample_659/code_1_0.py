def calculate_max_handstand_score():
    """
    Calculates the maximum possible score on MAG Floor with only handstand skills,
    based on the 2022 FIG Code of Points.
    """

    # Step 1: Identify all unique handstand press skills and their difficulty values (DV)
    # E=0.5, D=0.4, C=0.3, B=0.2, A=0.1
    skills = {
        "Manna press to handstand (E)": 0.5,
        "Nakayama press to handstand (E)": 0.5,
        "Straight body press through planche (D)": 0.4,
        "V-sit press to handstand (C)": 0.3,
        "Straddle sit press to handstand (C)": 0.3,
        "Hollow back press to handstand (B)": 0.2,
        "Press to handstand (A)": 0.1
    }

    # --- D-SCORE CALCULATION ---
    print("Calculating the Difficulty Score (D-score):\n")

    # Step 2: Sum the difficulty values of the available skills
    print("1. Sum of Skill Difficulty Values (DV):")
    difficulty_value_sum = 0
    for name, value in skills.items():
        print(f"   - {name}: {value:.1f}")
        difficulty_value_sum += value
    print(f"   ------------------------------------")
    print(f"   Total DV = {difficulty_value_sum:.1f}\n")

    # Step 3: Calculate Element Group (EG) bonus
    # Only EG I (non-acrobatic) is fulfilled. EG II, III, and IV (acrobatic groups) are not.
    element_group_bonus = 0.5
    print("2. Element Group (EG) Bonus:")
    print("   - EG I (Non-Acrobatic): Fulfilled (+0.5)")
    print("   - EG II, III, IV (Acrobatic): Not fulfilled (+0.0)")
    print(f"   Total EG Bonus = {element_group_bonus:.1f}\n")

    # Step 4: Calculate Connection Bonus (CB)
    # The max bonus for connecting non-acrobatic (EG I) elements is for a series
    # of 3 or more skills (one must be C-value or higher).
    connection_bonus = 0.2
    print("3. Connection Bonus (CB):")
    print("   - Connecting 3+ EG I elements gives a bonus of 0.2.")
    print(f"   Total CB = {connection_bonus:.1f}\n")

    # Step 5: Calculate total D-score
    d_score = difficulty_value_sum + element_group_bonus + connection_bonus
    print("---")
    print("Total D-Score = DV + EG Bonus + CB")
    print(f"Total D-Score = {difficulty_value_sum:.1f} + {element_group_bonus:.1f} + {connection_bonus:.1f} = {d_score:.1f}")
    print("---")

    # --- E-SCORE AND DEDUCTIONS ---
    
    # Step 6: Determine E-score
    e_score = 10.0
    print("\nCalculating the Execution Score (E-score):\n")
    print("Perfect execution is assumed.")
    print(f"E-score = {e_score:.1f}\n")

    # Step 7: Identify deductions
    # A 0.3 deduction applies if the last element is not a valid acrobatic dismount.
    deduction = 0.3
    print("Neutral Deductions:\n")
    print("A routine must end with an acrobatic dismount.")
    print(f"Deduction for no dismount = {deduction:.1f}\n")

    # Step 8: Calculate the final score
    final_score = d_score + e_score - deduction
    print("---")
    print("Final Score Calculation:\n")
    print("Final Score = (D-Score + E-Score) - Deductions")
    print(f"Final Score = ({d_score:.1f} + {e_score:.1f}) - {deduction:.1f}")
    print(f"Highest Possible Score = {final_score:.1f}")
    
    return final_score

# Run the calculation and store the final answer
final_answer = calculate_max_handstand_score()

# The final answer must be returned in the specified format
# print(f"\n<<<{final_answer:.1f}>>>")