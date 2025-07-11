import sympy

def solve_electostatics_problem():
    """
    Solves for the electric potential and field outside a conducting sphere
    in a uniform electric field using symbolic mathematics.
    """
    # Step 1 & 2: Define symbols and general form of potentials
    r, theta, R, E0, s1, s2 = sympy.symbols('r, theta, R, E_0, sigma_1, sigma_2')
    A, B = sympy.symbols('A, B')
    
    # P_1(cos(theta)) = cos(theta)
    cos_theta = sympy.cos(theta)
    
    # General potential solutions
    Phi_in = A * r * cos_theta
    Phi_out = -E0 * r * cos_theta + B * r**-2 * cos_theta
    
    # Step 3: Apply boundary conditions at r=R
    
    # Condition 3a: Continuity of potential
    # Phi_in(R, theta) = Phi_out(R, theta)
    eq1 = sympy.Eq(Phi_in.subs(r, R), Phi_out.subs(r, R))
    
    # Condition 3b: Continuity of normal current density
    # sigma_1 * d(Phi_in)/dr = sigma_2 * d(Phi_out)/dr at r=R
    dPhi_in_dr = sympy.diff(Phi_in, r)
    dPhi_out_dr = sympy.diff(Phi_out, r)
    eq2 = sympy.Eq(s1 * dPhi_in_dr.subs(r, R), s2 * dPhi_out_dr.subs(r, R))
    
    # Step 4: Solve for coefficients A and B
    solution = sympy.solve([eq1, eq2], (A, B))
    A_sol = solution[A]
    B_sol = solution[B]
    
    # Step 5: Substitute coefficients to find final potential and field
    
    # Final potential outside the sphere
    Phi_out_final = Phi_out.subs(B, B_sol)
    
    # Final electric field outside the sphere (E = -grad(Phi))
    E_out_r = -sympy.diff(Phi_out_final, r)
    E_out_theta = -sympy.diff(Phi_out_final, theta) / r
    
    # Simplify the expressions
    Phi_out_final_simp = sympy.simplify(Phi_out_final)
    E_out_r_simp = sympy.simplify(E_out_r)
    E_out_theta_simp = sympy.simplify(E_out_theta)

    # --- Output the results ---
    print("The problem asks for the electric potential and field in the region outside the sphere (r > R).")
    print("The following results are derived symbolically and match Option B.\n")

    # Reconstruct expressions for clean printing
    coeff_B_div_E0 = B_sol / E0
    
    print("Electric Potential Phi(r, theta) for r > R:")
    print(f"Phi(r, theta) = -{E0} * (r - ({sympy.printing.pretty(coeff_B_div_E0)}) / r**2) * cos(theta)")
    
    # Compare with Option B to be sure
    Phi_B = -E0 * (r - (s1 - s2) * R**3 / ((s1 + 2*s2) * r**2)) * cos_theta
    if sympy.simplify(Phi_out_final_simp - Phi_B) == 0:
        print("This matches the potential in Option B.\n")
    else:
        print("This does NOT match the potential in Option B.\n")
        
    print("Electric Field E(r, theta) for r > R:")
    
    # E_r component
    factor_r = sympy.simplify(E_out_r_simp / (E0 * cos_theta))
    print(f"E_r = {E0} * ({sympy.printing.pretty(factor_r)}) * cos(theta)")
    
    # E_theta component
    factor_theta = sympy.simplify(E_out_theta_simp / (-E0 * sympy.sin(theta)))
    print(f"E_theta = -{E0} * ({sympy.printing.pretty(factor_theta)}) * sin(theta)")

    # Compare with Option B to be sure
    E_B_r = E0 * (1 + 2*(s1 - s2)*R**3 / ((s1 + 2*s2) * r**3)) * cos_theta
    E_B_theta = -E0 * (1 - (s1 - s2)*R**3 / ((s1 + 2*s2) * r**3)) * sympy.sin(theta)
    
    if sympy.simplify(E_out_r_simp - E_B_r) == 0 and sympy.simplify(E_out_theta_simp - E_B_theta) == 0:
        print("This matches the electric field in Option B.")
    else:
        print("This does NOT match the electric field in Option B.")


if __name__ == "__main__":
    solve_electostatics_problem()
