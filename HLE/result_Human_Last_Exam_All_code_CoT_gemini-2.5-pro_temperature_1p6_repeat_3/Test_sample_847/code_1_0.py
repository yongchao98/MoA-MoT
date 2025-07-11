import numpy as np
import math

def solve_sorting_cost():
    """
    Calculates the coefficients for the asymptotic minimal cost to sort an array.

    The problem is solved by finding the most cost-effective information-gathering
    strategy. There are two types of questions Bob can ask. We analyze the cost
    per bit of information for each.

    Type 1 (Comparison): Cost=2 for 1 bit of information. Cost/bit = 2.

    Type 2 (General Yes/No): Cost is 1 for "yes" and 3 for "no". We can optimize
    the questions to minimize the worst-case cost. This leads to a cost model
    C(K) = min_k(max(1 + C(k), 3 + C(K-k))). For large K, assuming
    C(K) = alpha * log2(K), we find an optimal splitting ratio r that
    solves the cubic equation r^3 + 2r^2 + r - 1 = 0. The optimal
    cost per bit, alpha, can be calculated from r.

    This script solves the cubic equation for r, then computes the coefficients
    of the asymptotic cost formula: Cost = alpha * n*log2(n) - beta * n.
    """

    # Coefficients of the polynomial r^3 + 2*r^2 + r - 1 = 0
    coeffs = [1, 2, 1, -1]

    # Find the roots of the polynomial
    roots = np.roots(coeffs)

    # Find the single positive real root
    r = 0.0
    for root in roots:
        # Check if the root is real and positive
        if root.imag == 0 and root.real > 0:
            r = root.real
            break

    # Calculate alpha, the cost per bit of information using the optimal strategy
    # alpha = 1 / log2(r + 1)
    alpha = 1 / math.log2(r + 1)

    # Compare alpha with the cost/bit of Type 1 questions (which is 2).
    # Since alpha < 2, using Type 2 questions is always better.

    # The total cost is alpha * log2(n!)
    # Using Stirling's approximation: log2(n!) ~ n*log2(n) - n*log2(e)
    # The cost formula is Cost = alpha * n*log2(n) - alpha*log2(e) * n
    # Let beta = alpha * log2(e)
    log2_e = 1 / math.log(2)
    beta = alpha * log2_e
    
    # Print the coefficients and the final equation.
    print(f"The positive real root of the characteristic equation is r = {r:.5f}")
    print(f"The optimal cost per bit of information is alpha = {alpha:.3f}")
    print(f"The coefficient for the second term is beta = alpha * log2(e) = {beta:.3f}")
    print("\nThe minimal cost for large n is given by the asymptotic formula:")
    print(f"Cost = {alpha:.3f} * n*log2(n) - {beta:.3f} * n")

solve_sorting_cost()

# The question asks for "the answer" up to 3 decimal places.
# This most likely refers to the primary coefficient, alpha, which determines
# the leading term of the cost.
final_answer = 1 / math.log2(np.roots([1, 2, 1, -1])[np.isreal(np.roots([1, 2, 1, -1])) & (np.roots([1, 2, 1, -1]) > 0)][0].real + 1)
print(f"\n<<<{final_answer:.3f}>>>")