import math

def calculate_critical_exponent_nu(dimension):
    """
    Calculates and explains the critical exponent ν for a G₄ (φ⁴) theory
    in a given spatial dimension.

    The value of ν depends critically on the dimension 'd'.
    """

    print(f"Analysis for spatial dimension d = {dimension}:")
    print("---------------------------------------------")

    if dimension > 4:
        # For dimensions above the upper critical dimension (d_c = 4),
        # mean-field theory gives the exact result.
        nu = 0.5
        print("  - Physical Regime: Above the upper critical dimension (d > 4).")
        print("  - Applicable Theory: Mean-Field Theory.")
        print("  - In this regime, the critical exponent ν is exactly 1/2.")
        print(f"  - Final Equation: ν = 1 / 2")
        print(f"  - Precise Value: ν = {nu}")

    elif dimension == 4:
        # At the upper critical dimension, mean-field theory is the baseline,
        # but with important logarithmic corrections to scaling behavior.
        nu = 0.5
        print("  - Physical Regime: The upper critical dimension (d = 4).")
        print("  - Applicable Theory: Renormalization Group.")
        print("  - The exponent ν is 1/2, but with multiplicative logarithmic corrections.")
        print(f"  - Base Value: ν = {nu}")

    elif dimension == 3:
        # For d=3, the system belongs to the 3D Ising universality class.
        # The value is known with high precision from numerical methods
        # like Monte Carlo simulations and conformal bootstrap.
        nu = 0.6301
        print("  - Physical Regime: Below the upper critical dimension (d = 3).")
        print("  - Applicable Theory: 3D Ising Universality Class.")
        print("  - This value is determined from high-precision numerical and experimental results.")
        print(f"  - Approximate Value: ν ≈ {nu}")

    elif dimension == 2:
        # For d=2, the theory maps to the 2D Ising model,
        # which was solved exactly by Lars Onsager.
        nu = 1.0
        print("  - Physical Regime: Below the upper critical dimension (d = 2).")
        print("  - Applicable Theory: Exact solution of the 2D Ising Model.")
        print("  - The critical exponent ν is known to be exactly 1.")
        print(f"  - Precise Value: ν = {nu}")

    elif 1 < dimension < 4:
        # For d slightly below 4, the epsilon expansion (expansion in ε = 4-d)
        # from the Renormalization Group gives an excellent approximation.
        epsilon = 4 - dimension
        # One-loop result: ν = 1/2 + ε/12
        term1 = 0.5
        term2 = epsilon / 12.0
        nu = term1 + term2
        print(f"  - Physical Regime: Below the upper critical dimension (d < 4).")
        print("  - Applicable Theory: Renormalization Group (Epsilon Expansion).")
        print(f"  - The first-order epsilon expansion (ε = 4 - d = {epsilon:.2f}) gives:")
        print(f"  - Final Equation: ν = 1/2 + ε/12")
        # Outputting each number in the final equation as requested
        print(f"  - Calculation: ν = {term1} + {epsilon:.4f}/12.0 = {term1} + {term2:.4f}")
        print(f"  - Approximate Value: ν ≈ {nu:.4f}")
        
    else:
        # For d<=1, continuous phase transitions of this type do not occur for
        # systems with short-range interactions.
        print(f"  - For d={dimension}, a finite-temperature phase transition is not expected.")
        print("  - The exponent ν is not relevant in this context.")

    print("\n")


if __name__ == '__main__':
    # Demonstrate the calculation for a range of key dimensions
    dimensions_to_test = [5, 4, 3.5, 3, 2, 1]
    for d in dimensions_to_test:
        calculate_critical_exponent_nu(d)
