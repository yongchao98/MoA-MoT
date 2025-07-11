def calculate_max_handstand_score():
    """
    Calculates the maximum possible score for a MAG floor routine
    containing only handstand skills, assuming perfect execution and
    available connections.
    """
    
    # 1. Select the top 10 unique handstand-related skills from the 2022 FIG CoP
    # Skills are (Name, Difficulty Value)
    skills = [
        ("Handstand with 2/1 turn", 0.5), # E-value
        ("Press to Japanese handstand", 0.4), # D-value
        ("Handstand with 3/2 turn", 0.4), # D-value
        ("Japanese handstand to L-sit", 0.3), # C-value
        ("Handstand with 1/1 turn (Healy)", 0.3), # C-value
        ("Press to handstand (straight arms/body)", 0.2), # B-value
        ("From handstand, lower to support scale", 0.2), # B-value
        ("Felge to handstand", 0.2), # B-value
        ("Handstand with 1/2 turn", 0.2), # B-value
        ("Handstand with holding 2s", 0.1) # A-value
    ]

    # 2. Calculate D-Score components
    
    # Difficulty Value (DV): Sum of the values of the 10 skills
    difficulty_value = sum(skill[1] for skill in skills)
    
    # Difficulty Group Requirements (DGR): All skills are in Group I
    dgr_value = 0.5
    
    # Connection Value (CV): Assumed available as per the prompt.
    # The highest value is for connecting a D+ hold to a C+ hold from Group I.
    connection_value = 0.3
    
    # Total D-Score
    d_score = difficulty_value + dgr_value + connection_value
    
    # 3. Define E-Score (perfect execution)
    e_score = 10.0
    
    # 4. Calculate Final Score
    final_score = d_score + e_score

    # 5. Print the breakdown and final equation
    print("To achieve the maximum score with only handstand skills, the gymnast's score is calculated as follows:\n")
    print("--- D-Score Calculation ---")
    print("1. Difficulty Value (DV) from 10 skills:")
    skill_values_str = " + ".join([f"{val:.1f}" for _, val in skills])
    print(f"   DV = {skill_values_str} = {difficulty_value:.1f}\n")

    print("2. Difficulty Group Requirements (DGR):")
    print(f"   Fulfilling only Group I (Non-acrobatic) = {dgr_value:.1f}\n")

    print("3. Connection Value (CV):")
    print(f"   Connecting two Group I hold elements = {connection_value:.1f}\n")
    
    print(f"Total D-Score = {difficulty_value:.1f} (DV) + {dgr_value:.1f} (DGR) + {connection_value:.1f} (CV) = {d_score:.1f}\n")
    
    print("--- E-Score ---")
    print(f"Perfect Execution Score = {e_score:.1f}\n")

    print("--- Final Score Calculation ---")
    print("Final Score = D-Score + E-Score\n")
    
    # Final equation format: D-Score components shown explicitly
    final_equation_str = (
        f"{final_score:.1f} = "
        f"({skill_values_str}) "  # Skill values
        f"+ {dgr_value:.1f} "     # DGR
        f"+ {connection_value:.1f} " # CV
        f"+ {e_score:.1f}"       # E-Score
    )
    print("Final Equation:")
    print(final_equation_str)

calculate_max_handstand_score()