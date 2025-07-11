import argparse
import sys

def calculate_fn_p(n, p):
    """
    Calculates the value of the polynomial f_n(p) using the derived formula
    and prints the calculation steps.

    The formula is f_n(p) = (p^n - (1-p)^n) / (2p - 1).

    Args:
        n (int): The index 'n' of the polynomial.
        p (float): The value 'p' to evaluate the polynomial at.
    """
    
    # Handle the special case p = 0.5, where the denominator is zero.
    # In this case, we use L'Hopital's rule on the formula, which gives n * p^(n-1).
    if p == 0.5:
        # The value is n * (0.5)^(n-1)
        result = n * (0.5**(n - 1))
        print(f"For the special case p=0.5, the formula simplifies to f_n(0.5) = n * 0.5^(n-1)")
        print(f"f_{n}(0.5) = {n} * 0.5^({n-1})")
        print(f"f_{n}(0.5) = {n} * {0.5**(n-1)}")
        print(f"f_{n}(0.5) = {result}")
        return

    # Calculate the components of the formula
    p_to_the_n = p**n
    one_minus_p = 1 - p
    one_minus_p_to_the_n = one_minus_p**n
    numerator = p_to_the_n - one_minus_p_to_the_n
    denominator = 2 * p - 1

    # This check is technically redundant if p is a float due to the p==0.5 check,
    # but it is good practice to prevent division by zero.
    if denominator == 0:
        print("Error: Division by zero. This case should have been handled.", file=sys.stderr)
        return

    result = numerator / denominator

    # Print the calculation steps, showing each number in the equation.
    print(f"The general formula is: f_n(p) = (p^n - (1-p)^n) / (2*p - 1)")
    print(f"\nSubstituting n={n} and p={p}:")
    
    # Step 1: Initial substitution
    print(f"f_{n}({p}) = ({p}^{n} - (1-{p})^{n}) / (2*{p} - 1)")
    
    # Step 2: Evaluate the terms in parentheses and the denominator expression
    print(f"f_{n}({p}) = ({p_to_the_n} - ({one_minus_p})^{n}) / ({denominator})")

    # Step 3: Evaluate the second power term
    print(f"f_{n}({p}) = ({p_to_the_n} - {one_minus_p_to_the_n}) / {denominator}")

    # Step 4: Evaluate the numerator
    print(f"f_{n}({p}) = {numerator} / {denominator}")
    
    # Step 5: Final result
    print(f"f_{n}({p}) = {result}")

def main():
    """
    Parses command-line arguments for n and p and calls the calculation function.
    """
    parser = argparse.ArgumentParser(
        description="Calculate f_n(p) for a given sequence of polynomials.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Example usage:\npython solution.py 5 3.0\npython solution.py 10 0.5"
    )
    parser.add_argument('n', type=int, help="The index 'n' of the polynomial (an integer).")
    parser.add_argument('p', type=float, help="The value 'p' to evaluate the polynomial at (a float).")
    
    # Show help if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
        
    args = parser.parse_args()
    
    calculate_fn_p(args.n, args.p)

if __name__ == "__main__":
    main()