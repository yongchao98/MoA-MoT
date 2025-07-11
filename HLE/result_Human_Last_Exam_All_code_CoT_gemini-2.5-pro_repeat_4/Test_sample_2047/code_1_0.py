import math
import argparse
import sys

def solve_l(d, lam):
    """
    Calculates the value of l(d, lambda) based on the derived formula.

    Args:
        d (int): The dimension, must be an integer >= 4.
        lam (float): The lambda parameter, must be a float >= 1.
    """
    # Validate inputs based on the problem definition
    if not isinstance(d, int) or d < 4:
        print("Error: d must be an integer greater than or equal to 4.", file=sys.stderr)
        return
    if not isinstance(lam, (int, float)) or lam < 1:
        print("Error: lambda must be a number greater than or equal to 1.", file=sys.stderr)
        return

    try:
        # Calculate intermediate dot products
        x1_dot_p0 = math.sqrt(3 / d)
        x2_dot_p0 = math.sqrt(2 / d)

        # Calculate theta1 and theta2
        theta1 = math.acos(x1_dot_p0)
        theta2 = math.acos(x2_dot_p0)
    except ValueError as e:
        print(f"Error during calculation, likely due to invalid d: {e}", file=sys.stderr)
        return

    # Calculate the three terms of the l(d, lambda) formula
    # Term 1: The part dependent on lambda
    term_lambda = (theta2**2 - theta1**2) / (2 * lam)

    # Term 2: The log of the ratio of thetas
    term_log_theta = math.log(theta1 / theta2)

    # Term 3: The log of the ratio of sines
    # 0.5 * log((1-2/d)/(1-3/d)) = 0.5 * log((d-2)/(d-3))
    term_log_sin = 0.5 * math.log((d - 2) / (d - 3))
    
    # The final result is the sum of the three terms
    result = term_lambda + term_log_theta + term_log_sin

    # Print the detailed calculation steps as requested
    print(f"Calculating l(d, λ) for d={d} and λ={lam}")
    print(f"Formula: l(d, λ) = (θ₂² - θ₁²)/(2λ) + ln(θ₁/θ₂) + 0.5 * ln((d-2)/(d-3))")
    print("-" * 30)

    print(f"Step 1: Calculate θ₁ and θ₂")
    print(f"  θ₁ = arccos(sqrt(3/{d})) = arccos({x1_dot_p0:.6f}) = {theta1:.6f} rad")
    print(f"  θ₂ = arccos(sqrt(2/{d})) = arccos({x2_dot_p0:.6f}) = {theta2:.6f} rad")
    print("-" * 30)

    print(f"Step 2: Calculate each term of the equation")
    print(f"  Term 1 (λ-dependent): ({theta2**2:.6f} - {theta1**2:.6f}) / (2 * {lam}) = {term_lambda:.6f}")
    print(f"  Term 2 (θ-ratio): ln({theta1:.6f} / {theta2:.6f}) = ln({theta1/theta2:.6f}) = {term_log_theta:.6f}")
    print(f"  Term 3 (sin-ratio): 0.5 * ln(({d}-2)/({d}-3)) = 0.5 * ln({(d-2)/(d-3):.6f}) = {term_log_sin:.6f}")
    print("-" * 30)

    print(f"Step 3: Sum the terms to get the final result")
    print(f"  l({d}, {lam}) = {term_lambda:.6f} + {term_log_theta:.6f} + {term_log_sin:.6f}")
    
    # Print the final answer in the specified format
    print(f"<<< {result} >>>")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""Calculates l(d, lambda) for the given problem.
                       Example: python your_script_name.py 5 2.0""",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("d", type=int, help="Dimension d (integer, d >= 4)")
    parser.add_argument("lam", type=float, help="Lambda parameter λ (float, λ >= 1)")

    args = parser.parse_args()
    solve_l(args.d, args.lam)