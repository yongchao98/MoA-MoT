import sys

def calculate_max_handstand_score():
    """
    Calculates the maximum possible score on MAG Floor with only handstand skills
    under the 2022 FIG Code of Points, assuming perfect execution.
    """
    
    # 1. Execution Score (E-Score)
    # Assumed to be perfect.
    e_score = 10.0
    
    # 2. Difficulty Score (D-Score) Calculation
    
    # a) Difficulty Value (DV) from the top 8 skills
    # Values are taken from the FIG CoP for non-acrobatic elements (Group I).
    # G=0.7, F=0.6, E=0.5, D=0.4, C=0.3
    top_8_skills_values = [0.7, 0.6, 0.5, 0.4, 0.4, 0.4, 0.3, 0.3]
    dv_score = sum(top_8_skills_values)
    
    # b) Element Group (EG) Requirements
    # Only Group I (non-acrobatic) is fulfilled.
    # Groups II (acro fwd), III (acro bwd), and IV (dismount) are not fulfilled.
    eg_values = [0.5, 0.0, 0.0, 0.0]
    eg_score = sum(eg_values)
    
    # c) Connection Value (CV)
    # Maximum 4 connections. Connections of two different non-acrobatic elements.
    # G(0.7)+F(0.6) -> D+D or higher = 0.2 CV
    # E(0.5)+D(0.4) -> D+D or higher = 0.2 CV
    # D(0.4)+D(0.4) -> D+D or higher = 0.2 CV
    # C(0.3)+C(0.3) -> C+C          = 0.1 CV
    cv_values = [0.2, 0.2, 0.2, 0.1]
    cv_score = sum(cv_values)
    
    # Total D-Score
    d_score = dv_score + eg_score + cv_score
    
    # 3. Final Score
    final_score = d_score + e_score
    
    # --- Print the detailed calculation ---
    print("Calculating the maximum score for a MAG Floor routine with only handstands:\n")
    
    print("--- D-Score Calculation ---")
    dv_str = " + ".join(map(str, top_8_skills_values))
    print(f"1. Difficulty Value (DV) = {dv_str} = {dv_score:.1f}")
    
    eg_str = " + ".join(map(str, eg_values))
    print(f"2. Element Group (EG) Credit = {eg_str} = {eg_score:.1f}")
    
    cv_str = " + ".join(map(str, cv_values))
    print(f"3. Connection Value (CV) = {cv_str} = {cv_score:.1f}")
    
    print("-" * 27)
    print(f"Total D-Score = {dv_score:.1f} (DV) + {eg_score:.1f} (EG) + {cv_score:.1f} (CV) = {d_score:.1f}\n")

    print("--- Final Score Calculation ---")
    print(f"Execution Score (E-Score) = {e_score:.1f}")
    print(f"Difficulty Score (D-Score) = {d_score:.1f}\n")
    print(f"Final Score = {d_score:.1f} (D-Score) + {e_score:.1f} (E-Score) = {final_score:.1f}")

    # For the final answer extraction
    sys.stdout.write(f"\n<<<{final_score:.1f}>>>")

calculate_max_handstand_score()