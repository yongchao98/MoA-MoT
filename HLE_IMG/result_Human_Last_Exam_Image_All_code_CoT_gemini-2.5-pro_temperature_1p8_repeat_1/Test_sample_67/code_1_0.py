import sympy as sp

def solve_electron_energy_problem():
    """
    Solves for the minimum energy of electron 1 for the described process by
    applying conservation laws and finding the threshold condition.
    """
    # --- 1. Define Symbols and Equations ---
    # Define symbolic variables. C represents the constant ħ²/(2m*).
    # E_g is the band gap energy. k are the wave vectors.
    E_g, C, k1_i, k2_i, k1_f, k2_f = sp.symbols('E_g C k1_i k2_i k1_f k2_f', real=True, positive=True)

    # --- 2. State Conservation Laws and E-k relations ---
    print("Step 1: State the conservation laws and E(k) relations.")
    print("---------------------------------------------------------")
    print("Conservation of Energy: E1_i + E2_i = E1_f + E2_f")
    print("Conservation of Momentum: k1_i + k2_i = k1_f + k2_f\n")
    print("The energies are given by:")
    print(f"E1_i = E_g + C*k1_i**2")
    print(f"E2_i = -C*k2_i**2")
    print(f"E1_f = E_g + C*k1_f**2")
    print(f"E2_f = E_g + C*k2_f**2\n")

    # --- 3. Formulate the Threshold Condition ---
    print("Step 2: Find the conditions for the minimum energy threshold.")
    print("---------------------------------------------------------")
    print("To find the minimum E1_i, we consider the most efficient process:")
    print("1. The final state energy is minimized when k1_f = k2_f.")
    print("2. The initial state of electron 2 (k2_i) can be any value.\n")
    print("From momentum conservation, k1_f = k2_f = (k1_i + k2_i) / 2.\n")
    
    # --- 4. Derive the relationship between initial wave vectors ---
    print("Step 3: Derive the relationship between k1_i and k2_i.")
    print("-------------------------------------------------------")
    print("Substituting the relations into the energy conservation equation yields:")
    # Derivation: C*k1_i**2 - C*k2_i**2 = E_g + 2*C*((k1_i + k2_i)/2)**2
    # Rearranging gives a quadratic equation for k2_i:
    # 3*k2_i**2 + (2*k1_i)*k2_i + (2*E_g/C - k1_i**2) = 0
    quadratic_in_k2i = sp.Eq(3*k2_i**2 + 2*k1_i*k2_i + (2*E_g/C - k1_i**2), 0)
    print("Quadratic equation for k2_i:", quadratic_in_k2i)
    
    print("\nFor a real solution for k2_i to exist, the discriminant must be >= 0.")
    
    # --- 5. Solve for the minimum k1_i ---
    # Discriminant: b^2 - 4ac >= 0
    a, b, c = 3, 2*k1_i, (2*E_g/C - k1_i**2)
    discriminant = b**2 - 4*a*c
    
    print(f"\nDiscriminant = (2*k1_i)^2 - 4 * (3) * (2*E_g/C - k1_i**2) = {sp.simplify(discriminant)}")
    print("Setting the discriminant >= 0 gives a condition on k1_i:")
    
    # Solve 16*k1_i**2 - 24*E_g/C >= 0 for k1_i**2
    inequality = sp.solve(sp.simplify(discriminant) >= 0, k1_i**2)
    print(f"k1_i^2 >= {inequality.rhs}")
    k1_i_sq_min = inequality.rhs

    # --- 6. Calculate the Minimum Energy ---
    print("\nStep 4: Calculate the minimum initial energy of electron 1.")
    print("----------------------------------------------------------")
    print("The minimum energy E1_i corresponds to the minimum k1_i^2.")
    print(f"E1_i_min = E_g + C * k1_i_min^2")
    
    E1_i_min = E_g + C * k1_i_sq_min
    final_result = sp.simplify(E1_i_min)
    numeric_coefficient = final_result.coeff(E_g)
    
    print("\n--- Final Answer ---")
    print(f"The minimum required energy for electron 1 is:")
    print(f"E1_min = E_g + C * ({k1_i_sq_min}) = {final_result}")
    
    val_p = float(numeric_coefficient.p)
    val_q = float(numeric_coefficient.q)
    
    print("\nFinal Equation:")
    print(f"E1_min = ({int(val_p)} / {int(val_q)}) * Eg = {val_p/val_q} * Eg")


# Execute the function to print the solution
solve_electron_energy_problem()