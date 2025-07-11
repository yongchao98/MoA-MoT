import math
from fractions import Fraction

def solve_titan_mass_problem():
    """
    Solves the Titan mass calculation problem by finding the optimal fractional
    approximation that minimizes error while adhering to 5-bit constraints.
    """

    # --- Setup ---
    # Define constants and initial values as fractions
    r = Fraction(1, 2)
    rho = Fraction(9, 10)
    four_thirds = Fraction(4, 3)

    # Calculate the true mass for error comparison
    pi_true = math.pi
    m_true = float(four_thirds * (r**3) * rho * pi_true)

    # Define valid 5-bit fractional approximations for pi
    pi_approximations = [
        Fraction(3, 1), Fraction(13, 4), Fraction(16, 5), Fraction(19, 6),
        Fraction(22, 7), Fraction(25, 8), Fraction(28, 9), Fraction(31, 10)
    ]

    # --- Calculation ---
    # Helper to check if a fraction's components fit in 5 bits (<= 31)
    def is_valid(f):
        return f.numerator <= 31 and f.denominator <= 31

    # Calculate the constant part of the equation: C = (4/3) * r^3 * rho
    r_cubed = r**3  # (1/2)^3 = 1/8
    # Step-by-step multiplication, which is valid as simplified results fit the constraint
    const_part = (four_thirds * r_cubed) * rho # (1/6) * (9/10) -> 3/20

    # Find the best pi approximation that leads to a valid result and minimal error
    best_result = {
        "mass_frac": None,
        "pi_frac": None,
        "error": float('inf')
    }

    for pi_frac in pi_approximations:
        # We only consider pi approximations that are themselves valid 5-bit fractions
        if not is_valid(pi_frac):
            continue
        
        mass_frac = const_part * pi_frac
        
        # Check if the final result is a valid 5-bit fraction
        if is_valid(mass_frac):
            error = abs(float(mass_frac) - m_true)
            if error < best_result["error"]:
                best_result = {
                    "mass_frac": mass_frac,
                    "pi_frac": pi_frac,
                    "error": error
                }

    # --- Output the Derivation and Answer ---
    final_pi = best_result['pi_frac']
    final_mass = best_result['mass_frac']
    
    print("Derivation of the mass of the rock:")
    print("Formula: Mass = 4/3 * π * r³ * ρ")
    print("\nStep 1: Substitute initial values and combine constants.")
    print(f"Mass = (4/3) * (1/2)³ * (9/10) * π")
    print(f"Mass = (4/3 * 1/8 * 9/10) * π")
    print(f"Mass = {const_part.numerator}/{const_part.denominator} * π")

    print("\nStep 2: Substitute the optimal 5-bit approximation for π.")
    print(f"The best valid approximation is π ≈ {final_pi.numerator}/{final_pi.denominator}.")
    
    print("\nStep 3: Final calculation.")
    # As requested, outputting each number in the final equation
    print(f"Mass = {const_part.numerator}/{const_part.denominator} * {final_pi.numerator}/{final_pi.denominator}")
    print(f"Mass = {final_mass.numerator}/{final_mass.denominator} kg")

    final_error = round(best_result['error'], 3)
    print(f"\nThe smallest absolute error is {final_error:.3f}.")
    print(f"<<<{final_error:.3f}>>>")

solve_titan_mass_problem()