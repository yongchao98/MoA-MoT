import numpy as np

def solve_integral_approximation():
    """
    This function develops and presents an analytical formula that approximates the integral
    I(eps) = integral from 0 to 15 of 1 / (eps + 9x^5 + 5x^6 + 9x^8) dx
    for the small epsilon regime.
    """

    # --- Explanation of the Method ---
    explanation = """
### Method Explanation

The integral $I(\epsilon) = \int_0^{15.00} \\frac{1}{\epsilon + 9.0 x^{5.0} + 5.0 x^{6.0} + 9.0 x^{8.0}} dx$ is analyzed in the small $\epsilon$ limit.

1.  **Singular Behavior**: As $\epsilon \to 0^+$, the term $\epsilon$ in the denominator vanishes. The function $f(x) = 9x^5 + 5x^6 + 9x^8$ is zero at $x=0$, causing the integrand to become very large near this point. This means the behavior of the integral is dominated by the region near $x=0$.

2.  **Asymptotic Expansion**: This type of problem can be solved using asymptotic expansions. The leading behavior is determined by the lowest power of $x$ in $f(x)$, which is $9x^5$. The balance between $\epsilon$ and $9x^5$ suggests a variable scaling $x \sim \epsilon^{1/5}$, which leads to an asymptotic series for $I(\epsilon)$ in powers of $\epsilon^{1/5}$.

3.  **Approximation Formula**: We approximate the integral with the first two terms of this series, which provides a highly accurate formula for small $\epsilon$:
    $I(\epsilon) \approx c_0 \cdot \epsilon^{p_0} + c_1 \cdot \epsilon^{p_1}$

4.  **Coefficient Calculation**: The coefficients $c_0, c_1$ and powers $p_0, p_1$ are derived from the expansion.
    - The powers are $p_0 = -4/5$ and $p_1 = -3/5$.
    - The coefficients are calculated from integrals of the form $\int_0^\infty \\frac{t^{s-1}}{1+t^n} dt = \\frac{\pi/n}{\sin(\pi s/n)}$.
      - $c_0 = \\frac{1}{9^{1/5}} \int_0^\infty \\frac{du}{1+u^5} = \\frac{1}{9^{1/5}} \\frac{\pi/5}{\sin(\pi/5)}$
      - $c_1 = -\\frac{2}{9^{7/5}} \int_0^\infty \\frac{u}{1+u^5} du = -\\frac{2}{9^{7/5}} \\frac{\pi/5}{\sin(2\pi/5)}$

The code below computes the numerical values for these coefficients.
"""
    print(explanation)

    # --- Numerical Calculation of Coefficients ---

    # Powers
    p0 = -4.0 / 5.0
    p1 = -3.0 / 5.0

    # Coefficient c0
    c0_term1 = 9.0**(-1.0/5.0)
    c0_term2 = (np.pi / 5.0) / np.sin(np.pi / 5.0)
    c0 = c0_term1 * c0_term2

    # Coefficient c1
    c1_term1 = -2.0 / (9.0**(7.0/5.0))
    c1_term2 = (np.pi / 5.0) / np.sin(2.0 * np.pi / 5.0)
    c1 = c1_term1 * c1_term2

    # --- Final Answer ---
    print("### Analytical Formula")
    print("The analytical formula that approximates I(epsilon) for small epsilon is:")
    print(f"I(eps) ≈ {c0:.6f} * eps^({p0}) + {c1:.6f} * eps^({p1})")
    print("\nWhere the numbers in the equation are:")
    print(f"First coefficient: {c0}")
    print(f"First power: {p0}")
    print(f"Second coefficient: {c1}")
    print(f"Second power: {p1}")


solve_integral_approximation()
<<<I(eps) ≈ 0.688643 * eps^(-0.8) - 0.060956 * eps^(-0.6)>>>