import sys
from fractions import Fraction

def compute_limit_value(k):
    """
    Computes the value of the limit for a given integer k >= 2.

    The problem is to compute lim_{m -> inf} ln(f(m))/ln(m), where f(m)
    is the guaranteed size of the largest K_{k,k}-free subgraph over any
    graph with m edges.

    This limit evaluates to (2k - 1) / (2k).
    """
    if not isinstance(k, int) or k < 2:
        raise ValueError("k must be an integer greater than or equal to 2.")
    
    numerator = 2 * k - 1
    denominator = 2 * k
    
    # Using Fraction for exact rational number representation
    result_fraction = Fraction(numerator, denominator)
    
    print(f"For k = {k}, the problem asks for the value of the limit alpha.")
    print("The formula for alpha is (2*k - 1) / (2*k).")
    print(f"Plugging in k = {k}:")
    print(f"alpha = (2*{k} - 1) / (2*{k})")
    print(f"alpha = ({numerator}) / ({denominator})")
    print(f"The value of the limit is: {result_fraction}")
    print(f"As a decimal: {float(result_fraction)}")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python your_script_name.py <k>")
        print("where <k> is an integer >= 2.")
        # Running with a default example value if no argument is provided
        print("\nRunning with a default example value k=3...")
        k_val = 3
        compute_limit_value(k_val)
    else:
        try:
            k_val = int(sys.argv[1])
            compute_limit_value(k_val)
        except ValueError as e:
            print(f"Error: {e}")
