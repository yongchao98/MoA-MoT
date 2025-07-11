def derive_field_equation():
    """
    This script outlines the derivation of the field equation for Symmetric
    Teleparallel Gravity with the action S = -c^4/(16piG) * integral(sqrt(-g)*Q) d^4x + S_m.
    It prints the final equation with all its numerical coefficients.
    """

    print("Step 1: Identify the theory and the corresponding f(Q) Lagrangian.")
    print("The gravitational action is S_g = -c^4/(16*pi*G) * integral(sqrt(-g) * Q d^4x).")
    print("In the standard f(Q) gravity formalism with action S = 1/(2*kappa^2) * integral(sqrt(-g)*f(Q)),")
    print(f"where kappa^2 = 8*pi*G/c^4, our theory corresponds to f(Q) = -Q.")
    print("-" * 50)

    print("Step 2: State the general field equation for f(Q) gravity.")
    print("The general equation is:")
    print("f'(Q) * ( (2/sqrt(-g)) * d_a(sqrt(-g)P^a_mn) + P_mab*Q_n^ab - 2*Q_abm*P^ab_n )")
    print("- (1/2)*f(Q)*g_mn = (8*pi*G/c^4) * T_mn")
    print("-" * 50)

    print("Step 3: Substitute f(Q) = -Q and f'(Q) = -1 into the equation.")
    print("f(Q) = -Q")
    print("f'(Q) = -1")
    print("Substituting these gives:")
    print("-1 * ( (2/sqrt(-g)) * d_a(sqrt(-g)P^a_mn) + P_mab*Q_n^ab - 2*Q_abm*P^ab_n )")
    print("- (1/2)*(-Q)*g_mn = (8*pi*G/c^4) * T_mn")
    print("-" * 50)

    print("Step 4: Simplify the resulting equation.")
    print("After distributing the signs, we get the final field equation.")
    print("\nThe terms of the final equation are:")

    term1 = "-2/sqrt(-g) * d_a(sqrt(-g) * P^a_mn)"
    term2 = "- P_mab * Q_n^ab"
    term3 = "+ 2 * Q_abm * P^ab_n"
    term4 = "+ (1/2) * Q * g_mn"
    rhs = "(8*pi*G/c^4) * T_mn"

    print(f"Term 1 (from variation of Q): {term1}")
    print(f"Term 2 (from variation of P): {term2}")
    print(f"Term 3 (from variation of P): {term3}")
    print(f"Term 4 (from variation of sqrt(-g)): {term4}")
    print(f"Right Hand Side (from matter action): {rhs}")
    print("-" * 50)
    
    print("Final assembled field equation:")
    # Note: Using mn for mu,nu and a,b for alpha,beta for ASCII clarity.
    final_equation = f"{term1} {term2} {term3} {term4} = {rhs}"
    print(final_equation)
    print("\nThis equation matches option E.")

if __name__ == '__main__':
    derive_field_equation()