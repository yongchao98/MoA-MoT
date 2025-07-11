import sympy

def generate_polynomial(n):
    """
    Generates and prints the polynomial f_n(p) using the derived formula.

    Args:
        n (int): The index of the polynomial in the sequence.
    """
    if not isinstance(n, int) or n <= 0:
        print("Please enter a positive integer for n.")
        return

    p = sympy.symbols('p')

    # The derived formula is f_n(p) = (p**n - (1-p)**n) / (2*p - 1)
    # The numbers in the formula are 1 (in 1-p), 2 and 1 (in 2p-1).
    # The 'n' is a variable.
    numerator = p**n - (1 - p)**n
    denominator = 2 * p - 1

    # Using sympy.simplify to perform the polynomial division and simplification
    f_n = sympy.simplify(numerator / denominator)

    print(f"The derived general formula is f_n(p) = (p**n - (1 - p)**n) / (2*p - 1)")
    print("-" * 20)
    print(f"For n = {n}, the polynomial f_{n}(p) is:")
    
    # Using sympy.expand to make sure it's in expanded polynomial form
    final_poly = sympy.expand(f_n)
    
    # Printing the result in a readable format
    sympy.pprint(final_poly)


if __name__ == "__main__":
    try:
        n_input = int(input("Enter an integer n > 0 to generate the polynomial f_n(p): "))
        generate_polynomial(n_input)
    except ValueError:
        print("Invalid input. Please enter an integer.")
