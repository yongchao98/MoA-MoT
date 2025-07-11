import sys

def solve_nim_probability(n, m):
    """
    Determines if the first player has a winning position with a probability > 50%
    in a 2D-Generalized NIM game on a random n x m binary matrix.

    This simplifies to checking if n + m > 2.

    The function prints the analysis and the result.
    """
    # The problem reduces to checking the inequality n + m > 2.
    # The final equation is n + m > 2. The numbers are n, m, and 2.
    final_equation_n = n
    final_equation_m = m
    final_equation_constant = 2

    print(f"For an n x m matrix of size {n} x {m}:")
    print("The first player has a winning position with probability > 50% if and only if the following condition is true:")
    print(f"{final_equation_n} + {final_equation_m} > {final_equation_constant}")

    # Perform the calculation
    sum_val = final_equation_n + final_equation_m
    is_winning_probable = sum_val > final_equation_constant
    
    print(f"Calculating the left side: {sum_val}")
    print(f"Checking the inequality: {sum_val} > {final_equation_constant}")

    if is_winning_probable:
        print("The condition is True.")
        result = 1
    else:
        print("The condition is False.")
        result = 0

    print(f"The value of f({n}, {m}) is: {result}")

if __name__ == '__main__':
    # You can change these values to test with different matrix sizes.
    # Default values are provided.
    if len(sys.argv) == 3:
        try:
            n_val = int(sys.argv[1])
            m_val = int(sys.argv[2])
        except ValueError:
            print("Invalid input. Please provide two integers for n and m.")
            sys.exit(1)
    else:
        # Example values if no command-line arguments are given
        n_val = 5
        m_val = 7
    
    solve_nim_probability(n_val, m_val)
    print("\n" + "="*20 + "\n")
    # Test the edge case where n+m=2
    solve_nim_probability(1,1)