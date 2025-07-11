import sympy
import sys

def generate_polynomial(n):
    """
    Generates the polynomial f_n(p) for a given integer n
    using the derived formula: f_n(p) = (p^n - (1-p)^n) / (2p-1).

    Args:
        n (int): The index of the polynomial in the sequence. Must be a positive integer.

    Returns:
        A sympy expression for the polynomial f_n(p).
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("The index n must be a positive integer.")

    p = sympy.Symbol('p')
    
    # The numerator p**n - (1-p)**n is always divisible by (2p-1).
    # sympy's simplify function handles the symbolic division correctly.
    # We must expand the numerator first for sympy to recognize the cancellation.
    numerator = sympy.expand(p**n - (1-p)**n)
    denominator = 2*p - 1
    
    fn_polynomial = sympy.simplify(numerator / denominator)
    return fn_polynomial

def main():
    """
    Main function to parse the command-line argument and print the result.
    """
    print("This script generates the polynomial f_n(p) using the formula:")
    print("f_n(p) = (p**n - (1-p)**n) / (2*p-1)\n")
    
    # Check for command-line argument for n
    if len(sys.argv) == 2:
        try:
            n_value = int(sys.argv[1])
        except ValueError:
            print(f"Error: Invalid input '{sys.argv[1]}'. Please provide a positive integer for n.", file=sys.stderr)
            sys.exit(1)
    else:
        # If no argument is given, demonstrate with an example value.
        print("Usage: python your_script_name.py <n>")
        print("No value for n provided. Demonstrating for n=9...\n")
        n_value = 9
    
    try:
        # Generate the polynomial
        polynomial = generate_polynomial(n_value)
        
        # Format the output string to clearly show the resulting polynomial
        # Sympy's default string representation is clean and shows all numbers.
        final_output = f"f_{n_value}(p) = {polynomial}"
        
        print(final_output)

    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()