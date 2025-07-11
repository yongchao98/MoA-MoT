import sympy as sp

def solve_absorption_equation():
    """
    This function formulates and prints the equations for the absorption
    cross-section for a molecular chain under two conditions:
    a) No intermolecular interaction.
    b) Nearest-neighbor interaction (Frenkel exciton model).
    """

    # Define the symbolic variables for the physical quantities
    # Using sp.Symbol allows us to build and display mathematical equations
    omega = sp.Symbol('omega', real=True, positive=True)      # Laser frequency
    omega_0 = sp.Symbol('omega_0', real=True, positive=True)  # Molecular transition frequency
    mu_0_abs = sp.Symbol('|μ₀|', real=True, positive=True) # Magnitude of molecular transition dipole moment
    tau = sp.Symbol('τ', real=True, positive=True)           # Gaussian pulse duration
    J = sp.Symbol('J', real=True)                           # Nearest-neighbor coupling energy
    hbar = sp.Symbol('ħ', real=True, positive=True)         # Reduced Planck's constant
    K = sp.Symbol('K', real=True, positive=True)            # A constant of proportionality

    # --- Case a) No interaction between molecules ---
    # The absorption cross-section per molecule is simply that of an isolated molecule.
    # The transition occurs at frequency ω₀ with a strength proportional to |μ₀|².
    
    sigma_a = K * mu_0_abs**2 * omega_0 * tau * sp.exp(-(omega - omega_0)**2 * tau**2)

    print("Case a) No interaction between molecules:")
    print("The equation for the absorption cross-section per molecule, σ_a(ω), is:")
    sp.pprint(sigma_a, use_unicode=True)
    print("\nWhere:")
    print(f"  {omega} (omega): Laser frequency")
    print(f"  {omega_0} (omega_0): Intrinsic transition frequency of a single molecule")
    print(f"  {mu_0_abs} (|μ₀|): Magnitude of the molecular transition dipole moment")
    print(f"  {tau} (tau): Duration of the Gaussian laser pulse")
    print(f"  {K} (K): A constant of proportionality incorporating fundamental constants")
    print("-" * 60)


    # --- Case b) Nearest-neighbor interaction ---
    # The interaction creates delocalized exciton states. Due to a selection rule,
    # only the k=0 exciton is optically active (it can absorb light).
    
    # The transition frequency is shifted by the coupling energy J.
    omega_exc = omega_0 + 2 * J / hbar

    # The absorption cross-section per molecule has its peak shifted to the new
    # exciton frequency. The oscillator strength is conserved and concentrated here.
    
    sigma_b = K * mu_0_abs**2 * omega_exc * tau * sp.exp(-(omega - omega_exc)**2 * tau**2)

    print("Case b) Nearest-neighbor interaction (Frenkel Exciton Model):")
    print("The equation for the absorption cross-section per molecule, σ_b(ω), is:")
    sp.pprint(sigma_b, use_unicode=True)

    print("\nWhere the symbols are as defined above, with the addition of:")
    print(f"  {J} (J): Nearest-neighbor coupling energy")
    print(f"  {hbar} (ħ): Reduced Planck's constant")
    print("\nThe key physical consequences of the interaction are:")
    print("1. The absorption peak is shifted from ω₀ to the new exciton frequency, ω_exc:")
    print("   ω_exc = ")
    sp.pprint(omega_exc, use_unicode=True)
    print("2. The entire absorption strength is concentrated at this single, collective transition.")
    print("-" * 60)
    
    # The prompt asks to output each number in the final equation. As the equation is
    # symbolic, we can break it down into its constituent mathematical parts.
    print("Dissecting the final equation for the interacting case (Case b):")
    print("The final expression is a product of several factors:")
    print("\n1. Proportionality Constant:")
    sp.pprint(K, use_unicode=True)
    print("\n2. Squared Molecular Transition Dipole Moment (determines intrinsic strength):")
    sp.pprint(mu_0_abs**2, use_unicode=True)
    print("\n3. Exciton Transition Frequency (determines peak location and affects intensity):")
    sp.pprint(omega_exc, use_unicode=True)
    print("\n4. Pulse Duration (affects peak intensity):")
    sp.pprint(tau, use_unicode=True)
    print("\n5. Gaussian Lineshape Factor (determines the shape and width of the peak):")
    sp.pprint(sp.exp(-(omega - omega_exc)**2 * tau**2), use_unicode=True)

if __name__ == '__main__':
    solve_absorption_equation()