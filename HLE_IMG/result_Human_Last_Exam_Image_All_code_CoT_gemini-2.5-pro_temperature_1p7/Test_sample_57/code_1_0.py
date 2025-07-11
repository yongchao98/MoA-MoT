def main():
    """
    This function provides the derived formula for the number of ways to tile the given shape.
    The formula is expressed in terms of Fibonacci numbers F_n.
    """
    
    # The derived expression for the number of ways, a_n, is 2 * F_n * F_{n-1}.
    
    # We will now print the components of this final expression as requested.
    
    # The constant number in the equation:
    constant_factor = 2
    
    # The symbolic parts of the equation:
    fibonacci_term_n = "F_n"
    fibonacci_term_n_minus_1 = "F_{n-1}"
    
    print("The number of ways to fill the shape is given by the expression:")
    print(f"{constant_factor} * {fibonacci_term_n} * {fibonacci_term_n_minus_1}")
    
    # As per the instructions, here is each numerical component of the final equation printed out.
    print("\nBreaking down the formula:")
    print(f"The constant number is: {constant_factor}")
    print(f"The first Fibonacci term is: {fibonacci_term_n}")
    print(f"The second Fibonacci term is: {fibonacci_term_n_minus_1}")

if __name__ == "__main__":
    main()
