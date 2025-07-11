import numpy as np
from scipy.integrate import quad

def analyze_borwein_problem():
    """
    This script computes the values of Borwein integrals and related quantities
    to help determine which of the statements in the problem are correct.

    Plan:
    1. Define the function to be integrated, which is a product of sinc functions.
    2. For n from 1 to 8, calculate the sum S_n = Sum_{k=2..n} (1/k) which determines if I_n = π/2.
    3. For each n, numerically calculate I_n using high-precision integration.
    4. Present the results in a table for clear analysis.
    5. Perform and print specific calculations needed to evaluate certain statements.
    """

    # Step 1: Define the function to be integrated.
    def my_sinc(x):
        """Computes sin(x)/x, handling the singularity at x=0."""
        if x == 0.0:
            return 1.0
        return np.sin(x) / x

    def integrand_product(x, n):
        """Computes the product term in the Borwein integral for a given n."""
        product = 1.0
        for k in range(1, n + 1):
            product *= my_sinc(x / k)
        return product

    # Step 2 & 3: Compute the integrals and the Borwein condition sum.
    pi_half = np.pi / 2
    sum_condition = 0.0
    integrals = []

    print("Numerical Evaluation of Borwein Integrals I_n")
    print("-" * 80)
    print(f"{'n':>3} | {'Sum_{k=2..n} 1/k':>15} | {'I_n (numerical)':>25} | {'I_n - π/2':>25}")
    print("-" * 80)

    for n in range(1, 9):
        # Calculate the sum for the condition check
        if n >= 2:
            sum_condition += 1.0 / n
        
        # Use scipy.integrate.quad for high-precision numerical integration
        integral_val, error = quad(integrand_product, 0, np.inf, args=(n,), epsabs=1e-14, epsrel=1e-14)
        integrals.append(integral_val)
        
        sum_str = f"{sum_condition:.5f}" if n > 1 else "0.00000"
        print(f"{n:3d} | {sum_str:>15} | {integral_val:25.20f} | {(integral_val - pi_half):25.20f}")
    print("-" * 80)

    # Step 5: Print analysis relevant to specific statements.
    # For Statement F: Check if |I_5 - π/2| < 10⁻⁵
    i_5 = integrals[4] # n=5 is at index 4
    abs_diff_n5 = abs(i_5 - pi_half)
    print("\nAnalysis for Statement F (n=5):")
    # Outputting the numbers in the final evaluation for F
    print(f"|I_5 - π/2| = |{i_5:.10f} - {pi_half:.10f}| = {abs_diff_n5:.10f}")
    is_less = abs_diff_n5 < 1e-5
    print(f"Is {abs_diff_n5:.10f} < 10⁻⁵ (0.00001)? This is {is_less}.")

    # For Statement J: Check if the integrand is always positive for n=4.
    n_j = 4
    # Test a point where sinc(x) is negative. x > pi
    x_test = 3.5 
    integrand_val_at_test = integrand_product(x_test, n_j)
    # Outputting the numbers in the final evaluation for J
    print("\nAnalysis for Statement J (n=4):")
    print(f"The integrand for n=4 at x = {x_test} is {integrand_val_at_test:.4f}, which is not positive.")

analyze_borwein_problem()