import sympy

def analyze_soliton_stability():
    """
    Analyzes the stability of a localized soliton in 3D using a scaling argument
    (Derrick's theorem) for a Hamiltonian with exchange and DMI terms.
    """
    # Define symbolic variables for the scaling parameter lambda,
    # the dimension d, and the energy components E_ex and E_DMI.
    lam, E_ex, E_DMI, d = sympy.symbols('lambda E_ex E_DMI d')

    # The problem is in 3D
    d_val = 3

    print("This script analyzes the stability of a localized soliton in 3D using a scaling argument.")
    print(f"The energy functional is E = integral[ A*(grad m)^2 + D*m.curl(m) ] dV.")
    print("Let's consider a localized static solution and test its stability by scaling its size by a factor lambda.\n")

    # Construct the expression for the total energy E(lambda) of a scaled configuration.
    # The exchange term (grad m)^2 scales as lambda^(-2). The volume element dV scales as lambda^d.
    # So the exchange energy E_ex scales as lambda^(d-2).
    # The DMI term m.curl(m) scales as lambda^(-1). The volume element dV scales as lambda^d.
    # So the DMI energy E_DMI scales as lambda^(d-1).
    E = E_ex * lam**(d-2) + E_DMI * lam**(d-1)
    
    print("The energy of the scaled configuration E(lambda) in d dimensions is:")
    print(f"E(lambda) = E_ex * lambda^({d-2}) + E_DMI * lambda^({d-1})")

    # Substitute d=3
    E_3d = E.subs(d, d_val)
    print(f"\nFor d = {d_val} (3D), the energy is:")
    print(f"E(lambda) = {E_3d}\n")

    # Calculate the first derivative of E(lambda) with respect to lambda.
    dE_dlam = sympy.diff(E, lam)
    
    # Evaluate the first derivative at lambda=1 and set it to zero for a stationary solution.
    stationary_condition_lhs = dE_dlam.subs(lam, 1)
    stationary_condition_3d_lhs = stationary_condition_lhs.subs(d, d_val)

    print("For a stationary solution, the first derivative must be zero at lambda = 1:")
    print(f"dE/dlambda |_(lambda=1, d={d_val}) = {stationary_condition_3d_lhs} = 0")
    
    # This equation gives the relationship between the energy terms for a stationary solution.
    # We solve it for E_DMI.
    stationary_eq = sympy.Eq(stationary_condition_3d_lhs, 0)
    relation = sympy.solve(stationary_eq, E_DMI)
    print(f"This implies the following relation must hold: E_DMI = {relation[0]}\n")

    # Calculate the second derivative of E(lambda) with respect to lambda.
    d2E_dlam2 = sympy.diff(dE_dlam, lam)
    
    # Evaluate at lambda=1. This determines stability.
    stability_expr = d2E_dlam2.subs(lam, 1)
    stability_expr_3d = stability_expr.subs(d, d_val)

    print("To check for stability, we examine the second derivative at lambda = 1:")
    print(f"d^2E/dlambda^2 |_(lambda=1, d={d_val}) = {stability_expr_3d}")
    print("For the soliton to be stable, this value must be non-negative (>= 0).\n")

    # Substitute the relation from the stationary condition into the stability expression.
    stability_final_expr = stability_expr_3d.subs(E_DMI, relation[0])
    
    print("Now, we substitute the relation E_DMI = -E_ex/2 into the stability condition:")
    print(f"Stability Check = {stability_expr_3d}")
    # The expression 2*E_DMI gets the substitution
    final_equation_lhs = 2 * relation[0]
    print(f"                = {final_equation_lhs}")

    # Final analysis and conclusion
    print("\nThe final equation for the stability check is:")
    print(f"d^2E/dlambda^2 |_(stationary) = {final_equation_lhs}")
    print("\nThe exchange energy, E_ex = integral[A*(grad m)^2]dV, must be positive (E_ex > 0) for any non-uniform soliton, as the stiffness A is positive.")
    print(f"Therefore, the second derivative is -1 * E_ex, which is strictly negative.")
    print("A negative second derivative indicates that the stationary point is an energy maximum with respect to scaling, not a minimum.")
    print("This means the soliton is unstable against collapse (shrinking) or expansion.\n")
    
    print("--- FINAL CONCLUSION ---")
    print("No, it is not possible to stabilize a localized soliton in 3D with only Heisenberg exchange and Dzyaloshinskii-Moriya interaction.")

if __name__ == '__main__':
    analyze_soliton_stability()
