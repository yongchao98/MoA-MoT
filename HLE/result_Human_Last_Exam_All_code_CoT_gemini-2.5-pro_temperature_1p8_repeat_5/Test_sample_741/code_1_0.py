import math
from scipy.optimize import brentq

def sum_series(v):
    """
    Computes the sum S_v = sum_{i=0 to inf} 1 / (Gamma(v+i+1) * i!)
    This corresponds to the user's formula with v = x-1.
    """
    total = 0.0
    # The loop needs to run for enough iterations for the sum to converge.
    # 150 is a safe upper limit as the terms decrease very rapidly.
    for i in range(150):
        term = 0.0
        try:
            arg = v + i + 1.0
            # Gamma is infinite at non-positive integers. The inverse term is 0.
            if arg <= 0 and abs(arg - round(arg)) < 1e-9:
                term = 0.0
            else:
                term = 1.0 / (math.factorial(i) * math.gamma(arg))
        except (ValueError, OverflowError):
            # Handles cases where math.gamma or math.factorial might fail.
            term = 0.0
        
        total += term
        # Check for convergence to stop early.
        if i > 20 and abs(total) > 1e-12 and abs(term/total) < 1e-15:
            break
    return total

def solve_for_x():
    """
    Finds the largest x for which the summation converges to 0.
    """
    # By testing, we find the largest root for v is in the interval [-3.1, -3.0].
    # sum_series(-3.0) is positive, sum_series(-3.1) is negative.
    try:
        v_root = brentq(sum_series, -3.1, -3.0)
    except (ImportError, ValueError) as e:
        print("Error: SciPy library might be missing or root not in bracket.")
        print("Please install it using: pip install scipy")
        print(f"Error details: {e}")
        return

    # The problem asks for x, and we solved for v = x-1.
    x_root = v_root + 1

    print(f"The largest x value for which the summation converges to 0 is: {x_root:.4f}")
    
    # As requested, printing the numbers in the final equation.
    print("\nThe summation equation with this value of x is:")
    print(f"  Sum from i=0 to inf of 1 / (({x_root:.4f} + i - 1)! * i!) = 0\n")

    print("This can be expanded as:")
    
    total = 0.0
    output_str = ""
    for i in range(5):
        try:
            arg = x_root + i - 1
            term_val = 1.0 / (math.factorial(i) * math.gamma(arg + 1))
            total += term_val
            output_str += f"({term_val:+.4f})"
            if i < 4:
                output_str += " + "
        except (ValueError, OverflowError):
            output_str += "(error)"

    print(f"  {output_str} + ...")
    print(f"The partial sum of the first 5 terms is: {total:.4f}")
    
    # Final answer for parsing.
    # Respond in format {-a.bbb} as requested.
    print(f"\nFinal Answer Formatted: {{{x_root:.3f}}}")


solve_for_x()