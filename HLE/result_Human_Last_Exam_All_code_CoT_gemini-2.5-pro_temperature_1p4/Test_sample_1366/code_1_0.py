import math

def solve_voa_questions():
    """
    Solves the three questions about the vertex operator algebra V(p).
    """

    # Part (a): Decomposition of V(p)
    # The level k = -2 + 1/p is an admissible level for the affine Lie algebra associated with sl_2.
    # The representation category for L_k(sl_2) at these levels is known to be non-semisimple.
    # This means that modules are not always completely reducible (i.e., they don't always decompose into a direct sum of simple modules).
    # Therefore, a direct sum decomposition of V(p) is generally not possible.
    # However, it is known that there exist indecomposable modules that have a filtration (a sequence of submodules)
    # where the successive quotients are simple. This represents a different, weaker form of decomposition.
    answer_a_part1 = "No"
    answer_a_part2 = "Yes"
    answer_a = f"{answer_a_part1}, {answer_a_part2}"

    # Part (b): Top-level dimension of L(p)_n
    # The prompt explicitly defines the module L(p)_n as having the top-level rho_n.
    # It also defines rho_n as the (n+1)-dimensional irreducible sl_2-module.
    # Therefore, the dimension of the top-level of L(p)_n is n+1.
    answer_b = "n+1"

    # Part (c): Minimal conformal weight for p=2
    # The conformal weight h_n of the module L(p)_n is given by the Sugawara construction formula:
    # h_n = C2(rho_n) / (2 * (k + g)), where g is the dual Coxeter number (g=2 for sl_2)
    # and C2(rho_n) is the eigenvalue of the quadratic Casimir operator on rho_n.
    # For the (n+1)-dimensional representation rho_n of sl_2, C2(rho_n) = n(n+2)/2.
    # The level k = -2 + 1/p, so k+2 = 1/p.
    # Plugging these in: h_n = (n(n+2)/2) / (2 * (1/p)) = p*n*(n+2) / 4.
    # For p = 2, this simplifies to h_n = 2*n*(n+2) / 4 = n*(n+2)/2.
    # The decomposition is a sum over n >= 0.
    # The function h_n = n*(n+2)/2 is monotonically increasing for n >= 0.
    # Thus, the minimum value occurs at the smallest value of n, which is n=0.
    # h_0 = 0*(0+2)/2 = 0.
    # The minimal conformal weight is 0.
    p_val = 2
    n_min = 0
    # Equation for minimal conformal weight: h_0 = (p * n * (n + 2)) / 4
    min_weight = (p_val * n_min * (n_min + 2)) / 4
    answer_c = int(min_weight)

    # Print the final combined answer in the required format.
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve_voa_questions()