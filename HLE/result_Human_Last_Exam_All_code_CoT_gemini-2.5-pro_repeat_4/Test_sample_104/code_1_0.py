import numpy as np
from scipy.stats import chi2, ncx2
from scipy.integrate import quad
from scipy.optimize import brentq

def solve_for_f_statistic():
    """
    Calculates the minimum F-statistic for a given relative bias tolerance and confidence level.
    """
    
    # --- Step 1: Find the non-centrality parameter (λ) for 10% bias ---
    
    # For TSLS with one instrument (L=1), the maximal relative bias is given by
    # the expectation E[V / (λ + V)], where V follows a chi-squared distribution with df=1.
    
    def integrand(x, lmbda):
        """The function to be integrated: (x / (λ + x)) * pdf_chi2(x, df=1)."""
        return (x / (lmbda + x)) * chi2.pdf(x, df=1)

    def expected_relative_bias(lmbda):
        """Calculates the expected relative bias for a given non-centrality parameter λ."""
        if lmbda <= 0:
            return np.inf
        # Integrate from 0 to infinity
        result, _ = quad(integrand, 0, np.inf, args=(lmbda,))
        return result

    def equation_to_solve(lmbda):
        """This is the function we want to find the root of: E[bias] - 0.10 = 0."""
        return expected_relative_bias(lmbda) - 0.10

    # Solve for λ_0.10, the value of λ where the maximal relative bias is exactly 0.10.
    # We use a numerical solver and provide a search bracket [0.1, 100].
    try:
        lambda_01 = brentq(equation_to_solve, 0.1, 100)
    except ValueError:
        print("Could not find a root for lambda in the given interval.")
        return

    # --- Step 2: Find the 95% critical value for the F-statistic ---

    # With one instrument, the F-statistic is distributed as a non-central chi-squared
    # with df=1 and non-centrality parameter λ.
    df = 1
    
    # We need to find the critical value F_crit such that P(F > F_crit | λ = λ_0.10) = 0.05.
    # This is the 95th percentile of the ncx2(df=1, nc=λ_0.10) distribution.
    critical_value = ncx2.ppf(0.95, df=df, nc=lambda_01)
    
    # Round to the nearest whole number for the final answer.
    final_answer = round(critical_value)

    # --- Output the results and the logic ---
    print(f"To have 95% confidence that the relative TSLS bias is less than 10% with one instrument:")
    print("-" * 75)
    print(f"Step 1: Find the non-centrality parameter (λ) where the maximal relative bias is exactly 10%.")
    print(f"This is the solution to the equation: E[V/(λ+V)] = 0.10, where V is a χ²(1) random variable.")
    print(f"The calculated non-centrality parameter is λ ≈ {lambda_01:.4f}.")
    print("-" * 75)
    print(f"Step 2: Find the critical F-statistic based on this λ.")
    print(f"The F-statistic follows a non-central χ² distribution with df=1 and nc=λ.")
    print(f"The critical value is the 95th percentile of this distribution.")
    print(f"F_critical = ppf(0.95, df=1, nc={lambda_01:.4f}) ≈ {critical_value:.4f}.")
    print("-" * 75)
    print(f"Conclusion: The minimum F-statistic required, rounded to the nearest whole number, is {final_answer}.")

solve_for_f_statistic()