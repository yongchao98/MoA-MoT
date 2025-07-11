import numpy as np
import math
from scipy.integrate import quad
from scipy.special import erf

def solve_projection_area():
    """
    Calculates the average area of a projection of the n-dimensional cross-polytope,
    where n = 2k + 1 is an odd dimension provided by the user.
    """
    try:
        k_str = input("Enter the value of k (where n=2k+1): ")
        k = int(k_str)
        if k < 0:
            print("k must be a non-negative integer.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    n = 2 * k + 1

    # Define the integrand for M_n
    def M_n_integrand(x):
        return 1 - (erf(x / np.sqrt(2)))**n

    # Calculate M_n = E[max(|Z_i|)] by numerical integration
    M_n, error_M_n = quad(M_n_integrand, 0, np.inf)

    # Calculate the coefficient C = (sqrt(2) * pi^k * Gamma(k+1/2)) / (2 * (k!)^2)
    # Using log-gamma for numerical stability
    log_gamma_k_half = math.lgamma(k + 0.5)
    log_factorial_k = math.lgamma(k + 1)
    
    # Calculate log of C
    log_C = 0.5 * np.log(2) + k * np.log(np.pi) + log_gamma_k_half - np.log(2) - 2 * log_factorial_k
    C = np.exp(log_C)

    # Calculate the average area
    average_area = C * M_n

    # Print the full final equation for the user
    # Note: printing the raw formula as a string
    k_val = k
    n_val = 2 * k_val + 1
    
    print(f"For k = {k_val}, the dimension is n = {n_val}.")
    print("\nThe equation for the average area (A) is:")
    
    # Display formula with values plugged in
    gamma_k_half_val = math.gamma(k_val + 0.5)
    k_fact_val = math.factorial(k_val)

    # To build the output string for the formula
    numerator_str = f"sqrt(2) * pi^{k_val} * Gamma({k_val}+0.5)"
    denominator_str = f"2 * ({k_val}!)^2"
    integral_str = f"Integral[1 - erf(x/sqrt(2))^{n_val}, {{x, 0, inf}}]"

    # Print a more readable equation form
    print(f"A = ({numerator_str}) / ({denominator_str}) * {integral_str}")
    
    # Show the value of each part of the formula
    print("\nEvaluating each part:")
    numerator_val = np.sqrt(2) * (np.pi**k_val) * gamma_k_half_val
    denominator_val = 2 * (k_fact_val**2)
    
    print(f"Value of the coefficient ({numerator_str}) / ({denominator_str}):")
    print(f"= ({np.sqrt(2):.4f} * {np.pi**k_val:.4f} * {gamma_k_half_val:.4f}) / (2 * ({k_fact_val})^2)")
    print(f"= {numerator_val:.4f} / {denominator_val:.4f} = {C:.4f}")

    print(f"\nValue of the integral {integral_str}: {M_n:.4f}")

    print("\nFinal Result:")
    print(f"Average Area = {C:.4f} * {M_n:.4f} = {average_area}")


solve_projection_area()