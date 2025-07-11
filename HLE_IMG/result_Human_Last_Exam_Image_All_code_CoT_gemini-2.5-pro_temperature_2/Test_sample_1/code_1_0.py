import math

def calculate_reaction_rate():
    """
    Identifies the rate-determining step and calculates the reaction rate constant
    based on the provided Gibbs Free Energy diagram.
    """
    # 1. Define constants and given values
    kB = 1.380649e-23  # Boltzmann constant in J/K
    h = 6.62607015e-34   # Planck constant in J·s
    R = 8.31446      # Ideal gas constant in J/(mol·K)
    T = 298          # Temperature in Kelvin from the graph axis

    # Gibbs Free Energies from the diagram (in kJ/mol)
    G_R = -55
    G_TS_add = 10
    G_P5 = -22
    G_TS_eli = 18

    # 2. Identify the rate-determining step by comparing activation barriers
    delta_G_add_barrier = G_TS_add - G_R
    delta_G_eli_barrier = G_TS_eli - G_P5

    print("Step 1: Identify the rate-determining step (RDS) and its activation energy.")
    print("The rate of a multi-step reaction is determined by the step with the highest activation energy barrier.")
    print(f"Barrier for addition step (R -> P5 via TS-add): ΔG‡_add = G(TS-add) - G(R) = {G_TS_add} kJ/mol - ({G_R} kJ/mol) = {delta_G_add_barrier} kJ/mol.")
    print(f"Barrier for elimination step (P5 -> P via TS-eli): ΔG‡_eli = G(TS-eli) - G(P5) = {G_TS_eli} kJ/mol - ({G_P5} kJ/mol) = {delta_G_eli_barrier} kJ/mol.")
    
    if delta_G_add_barrier > delta_G_eli_barrier:
        rds = "addition step (R -> P5)"
        delta_G_activated_kJ = delta_G_add_barrier
    else:
        rds = "elimination step (P5 -> P)"
        delta_G_activated_kJ = delta_G_eli_barrier

    delta_G_activated_J = delta_G_activated_kJ * 1000 # Convert kJ/mol to J/mol
    
    print(f"\nThe highest barrier is {delta_G_activated_kJ} kJ/mol. Therefore, the {rds} is rate-determining.")

    # 3. Calculate the rate constant using the Eyring equation
    print("\nStep 2: Calculate the reaction rate constant (k) using the Eyring equation.")
    print("The Eyring equation is: k = (kB * T / h) * exp(-ΔG‡ / (R * T))")
    print("Substituting the numerical values into the equation:")
    
    # Calculate each part of the equation
    pre_exponential_factor = (kB * T) / h
    exponent_val = -delta_G_activated_J / (R * T)
    k_in_seconds = pre_exponential_factor * math.exp(exponent_val)
    
    print(f"k = (({kB:.6e} J/K * {T} K) / {h:.6e} J·s) * exp(-{delta_G_activated_J} J/mol / ({R:.5f} J/(mol·K) * {T} K))")

    # 4. Convert to hours^-1 and format to two significant figures
    k_in_hours = k_in_seconds * 3600

    print(f"\nThe calculated rate constant in s^-1 is: k = {k_in_seconds:.2f} s⁻¹")
    print(f"\nStep 3: Convert the rate constant to units of hours⁻¹.")
    print(f"k (in hours⁻¹) = {k_in_seconds:.2f} s⁻¹ * 3600 s/hr")
    print(f"The final reaction rate constant is: {k_in_hours:.2g} hours⁻¹")

calculate_reaction_rate()
<<<9.1e+03>>>