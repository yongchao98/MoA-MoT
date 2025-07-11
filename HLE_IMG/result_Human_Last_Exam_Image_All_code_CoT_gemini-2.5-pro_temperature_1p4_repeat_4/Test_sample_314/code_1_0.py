def solve_graphene_band_structure():
    """
    Analyzes four graphene band structure simulations to determine which corresponds to
    specific tight-binding parameter conditions.
    """

    # Step 1 & 2: Qualitative Analysis of Each Simulation from the images.
    # Energies are estimated at the Gamma point (k=0 for 3D plots).

    # Simulation 1 (Line plot):
    # E_conduction(Gamma) approx -1 eV, E_valence(Gamma) approx -15 eV.
    # Note: Even though both energies are negative, which is unusual for a standard plot with E_fermi=0,
    # we can still analyze the properties.
    s1_sign = '+' # Valence band |E|=15 is much wider than conduction band |E|=1
    bandwidth1 = abs(-1 - (-15))  # 14 eV
    # Asymmetry ratio |E_v/E_c|. A larger ratio implies a larger |s|.
    asymmetry_ratio1 = abs(-15 / -1.0) # 15.0

    # Simulation 2 (3D plot):
    # E_conduction(Gamma) approx +2.5 eV, E_valence(Gamma) approx -10 eV.
    s2_sign = '+' # Valence band is wider
    bandwidth2 = abs(2.5 - (-10)) # 12.5 eV
    asymmetry_ratio2 = abs(-10 / 2.5) # 4.0

    # Simulation 3 (3D plot):
    # E_conduction(Gamma) approx +5 eV, E_valence(Gamma) approx -15 eV.
    s3_sign = '+' # Valence band is wider
    bandwidth3 = abs(5 - (-15)) # 20.0 eV
    asymmetry_ratio3 = abs(-15 / 5.0) # 3.0

    # Simulation 4 (3D plot):
    # E_conduction(Gamma) approx +15 eV, E_valence(Gamma) approx -5 eV.
    s4_sign = '-' # Conduction band is wider
    bandwidth4 = abs(15 - (-5)) # 20.0 eV
    # Asymmetry for s<0 is |E_c/E_v|
    asymmetry_ratio4 = abs(15 / -5.0) # 3.0

    # Step 3: Match simulations to the conditions.

    # Condition 3: unique sign(s)
    # Simulation 4 is the only one with a negative sign for s.
    sim_for_cond3 = 4

    # Condition 1: minimum t
    # Hopping parameter 't' is proportional to the bandwidth. We find the simulation with the smallest bandwidth.
    # Bandwidths: Sim1=14, Sim2=12.5, Sim3=20, Sim4=20.
    # Simulation 2 has the minimum bandwidth.
    sim_for_cond1 = 2

    # We now evaluate simulations 1 and 3 for the remaining conditions 2 and 4.
    # Both have s > 0.

    # Condition 2: minimum |s|
    # The magnitude of 's' is related to the asymmetry. We look for the smallest asymmetry ratio among {1, 3}.
    # Ratios: asymmetry_ratio1=15.0, asymmetry_ratio3=3.0.
    # Simulation 3 has a smaller asymmetry ratio than Simulation 1.
    sim_for_cond2 = 3

    # Condition 4: maximum s
    # This means the largest positive s. It corresponds to the greatest asymmetry among {1, 2, 3}.
    # Ratios: R1=15.0, R2=4.0, R3=3.0.
    # Simulation 1 has the largest asymmetry ratio.
    sim_for_cond1_check = sim_for_cond1
    sim_for_cond4 = 1

    # Final Answer Assembly
    # The final answer is the sequence of simulation indices for conditions 1, 2, 3, and 4.
    answer = f"{sim_for_cond1}{sim_for_cond2}{sim_for_cond3}{sim_for_cond4}"

    print("Analysis of Graphene Band Structures:")
    print("-" * 35)
    print(f"Condition 1: Minimum hopping parameter (t), corresponds to minimum bandwidth.")
    print(f"-> Simulation {sim_for_cond1} has the smallest bandwidth ({bandwidth2} eV).")
    print("-" * 35)
    print(f"Condition 2: Minimum overlap magnitude (|s|), corresponds to the least asymmetry (for s>0).")
    print(f"-> Simulation {sim_for_cond2} has the smallest asymmetry ratio ({asymmetry_ratio3}).")
    print("-" * 35)
    print(f"Condition 3: Unique overlap sign (sign(s)).")
    print(f"-> Simulation {sim_for_cond3} is the only one with s<0 (wider conduction band).")
    print("-" * 35)
    print(f"Condition 4: Maximum overlap magnitude (s), corresponds to the most asymmetry (for s>0).")
    print(f"-> Simulation {sim_for_cond4} has the largest asymmetry ratio ({asymmetry_ratio1}).")
    print("-" * 35)
    print(f"The simulation indices ordered by the condition met (1, 2, 3, 4) are: {answer}")
    
solve_graphene_band_structure()