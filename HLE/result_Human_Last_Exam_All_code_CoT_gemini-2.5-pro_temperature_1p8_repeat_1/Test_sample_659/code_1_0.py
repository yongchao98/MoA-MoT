def calculate_max_handstand_score():
    """
    Calculates the highest possible MAG floor score with only handstand skills.

    The calculation is based on the 2022-2024 FIG Code of Points, assuming
    perfect execution and availability of all handstand skills and connections.
    """

    # --- 1. E-Score (Execution) ---
    # Assuming perfect execution as per the prompt.
    e_score = 10.0
    print(f"E-Score (Execution): {e_score} (for perfect execution)\n")

    # --- 2. D-Score (Difficulty) ---
    print("Calculating the D-Score (Difficulty)...")

    # --- 2a. Skill Difficulty Values ---
    # The D-Score includes the top 10 elements. We will select the 10 highest-value
    # distinct handstand-related skills from the MAG Floor Table of Elements.
    # Difficulty values: B=0.2, C=0.3, D=0.4, E=0.5, F=0.6
    
    skills = {
        "Press from Manna to handstand (hold 2s)": 0.6,  # F Value
        "Press from Planche to handstand (hold 2s)": 0.5, # E Value
        "Handstand with 1080-degree turn": 0.5,            # E Value
        "Flair to handstand (hold 2s)": 0.4,              # D Value
        "Straight arm press from straddle sit to handstand (hold 2s)": 0.4, # D Value
        "Handstand with 720-degree turn": 0.4,             # D Value
        "Guczoghy (Straddle L to Handstand, lower to L)": 0.4, # D Value
        "Press from V-sit to handstand (hold 2s)": 0.3,    # C Value
        "Handstand with 360-degree turn": 0.3,             # C Value
        "Press from stand/pike to handstand (hold 2s)": 0.2 # B Value
    }

    skill_difficulty_sum = sum(skills.values())
    
    print(f"Top 10 Handstand Skills Sum: {skill_difficulty_sum:.1f}")
    # print the individual skills
    for name, val in skills.items():
        print(f"- {name}: {val:.1f}")
    
    print("")

    # --- 2b. Element Group Requirements (EGR) ---
    # There are 4 element groups. Bonus is 0.5 for each one fulfilled.
    # Group I: Non-acrobatic (fulfilled by handstands)
    # Group II: Acrobatic Forward (not fulfilled)
    # Group III: Acrobatic Backward (not fulfilled)
    # Group IV: Acrobatic Sideward (not fulfilled)
    egr_bonus = 0.5  # Only Group I is fulfilled.
    print(f"Element Group Requirement (EGR) Bonus: {egr_bonus:.1f} (for Group I only)")
    
    # --- 2c. Connection Value (CV) ---
    # A maximum of 3 connections can be awarded.
    # Connection of D+D (or higher) non-acrobatic elements gives a 0.2 bonus.
    # With F, E, E, D, D, D, D elements, three D+D connections are possible.
    cv_per_connection = 0.2
    num_connections = 3
    cv_bonus = cv_per_connection * num_connections
    print(f"Connection Value (CV) Bonus: {cv_bonus:.1f} (from 3 connections of D-level or higher elements)\n")

    # --- Total D-Score Calculation ---
    d_score = skill_difficulty_sum + egr_bonus + cv_bonus
    print(f"Total D-Score = {skill_difficulty_sum:.1f} (Skills) + {egr_bonus:.1f} (EGR) + {cv_bonus:.1f} (CV)")
    print(f"Total D-Score = {d_score:.1f}\n")
    
    # --- Final Score Calculation ---
    final_score = d_score + e_score
    
    print("--- Final Score Calculation ---")
    print(f"Final Score = {d_score:.1f} (D-Score) + {e_score:.1f} (E-Score)")
    print(f"Highest Possible Score = {final_score:.1f}")
    
    return final_score

# Run the calculation and get the final answer
highest_score = calculate_max_handstand_score()
print(f"\n<<< {highest_score:.1f} >>>")