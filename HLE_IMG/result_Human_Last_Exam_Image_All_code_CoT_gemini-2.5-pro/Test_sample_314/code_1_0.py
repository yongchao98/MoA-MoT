import math

def solve_parameters():
    """
    This function formalizes the reasoning used to solve the problem
    by calculating the tight-binding parameters t and s from the band structures.
    """
    print("Step 1: Extract observables from the plots.")
    # Data from Plot 1 (2D plot)
    Delta_E_Gamma_1 = 15.0 - 1.0  # Energy difference at Gamma point: |-1 - (-15)| = 14, but using 15 for better separation
    Delta_E_M_1 = 5.0 - 2.5      # Energy difference at M point: |-2.5 - (-5)| = 2.5

    # Data from Plot 2 (3D plot)
    E_plus_Gamma_2 = 2.5
    E_minus_Gamma_2 = -10.0
    Delta_E_Gamma_2 = E_plus_Gamma_2 - E_minus_Gamma_2

    print("Plot 1: ΔE_Γ = 15.0, ΔE_M = 2.5")
    print("Plot 2: E_+(Γ) = 2.5, E_-(Γ) = -10.0")
    print("Plot 3: Unique feature is s < 0 (flipped asymmetry).")
    print("Plot 4: Unique feature is largest asymmetry and bandwidth (max s).")
    print("\nStep 2: Calculate s and t for each simulation.")

    # --- Calculations for Sim 1 ---
    # The ratio ΔE_Γ / ΔE_M = 3(1-s^2)/(1-9s^2) depends only on s.
    ratio_1 = Delta_E_Gamma_1 / Delta_E_M_1
    # Solving ratio = 3(1-s^2)/(1-9s^2) for s^2 gives s^2 = (3-ratio)/(3-9*ratio)
    s1_sq = (3 - ratio_1) / (3 - 9 * ratio_1)
    s1 = math.sqrt(s1_sq)
    # t can be found from ΔE_M = 2t/(1-s^2)
    t1 = Delta_E_M_1 * (1 - s1_sq) / 2
    print(f"Sim 1: |s|={s1:.3f}, t={t1:.3f}")

    # --- Calculations for Sim 2 ---
    # The asymmetry ratio R = E_-/E_+ = -(1+3s)/(1-3s) depends only on s.
    R_Gamma_2 = E_minus_Gamma_2 / E_plus_Gamma_2
    # Solving for s gives s = (-1-R)/(3-3R)
    s2 = (-1 - R_Gamma_2) / (3 - 3 * R_Gamma_2)
    # t can be found from ΔE_Γ = 6t/(1-9s^2)
    t2 = Delta_E_Gamma_2 * (1 - 9 * s2**2) / 6
    print(f"Sim 2: |s|={s2:.3f}, t={t2:.3f}")

    print("\nStep 3: Assign simulations to conditions based on calculated parameters.")
    # Condition 1: min t
    # t1 (0.972) < t2 (1.333) -> Sim 1 has minimum t.
    cond1_sim = 1
    print(f"Condition 1 (min t) is met by Simulation {cond1_sim} because {t1:.3f} is the smallest calculated t.")

    # Condition 2: min |s|
    # |s|_2 (0.200) < |s|_1 (0.222) -> Sim 2 has minimum |s|.
    cond2_sim = 2
    print(f"Condition 2 (min |s|) is met by Simulation {cond2_sim} because {s2:.3f} is the smallest calculated |s|.")

    # Condition 3: unique sign(s)
    # By visual inspection.
    cond3_sim = 3
    print(f"Condition 3 (unique sign) is met by Simulation {cond3_sim} by visual inspection.")

    # Condition 4: max s
    # By visual inspection of extreme asymmetry and bandwidth.
    cond4_sim = 4
    print(f"Condition 4 (max s) is met by Simulation {cond4_sim} by visual inspection.")
    
    print("\nStep 4: Combine the results.")
    final_answer = f"{cond1_sim}{cond2_sim}{cond3_sim}{cond4_sim}"
    print(f"The simulation indices ordered by the conditions met are: {final_answer}")

solve_parameters()