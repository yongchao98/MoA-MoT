def display_q_derivative(n):
    """
    Calculates and prints the formula for the q-derivative of T^n for a non-negative integer n.
    The output includes all numerical components of the resulting equation.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: n must be a non-negative integer.")
        return

    # Handle the case for n=0.
    # The equation is nabla_q(T^0) = 0. The numbers involved are 0 and 0.
    if n == 0:
        print(f"nabla_q(T^{n}) = 0")
        return

    # Handle the case for n=1.
    # The equation is nabla_q(T^1) = 1. The numbers involved are 1 and 1.
    if n == 1:
        print(f"nabla_q(T^{n}) = 1")
        return

    # Handle the general case for n > 1.
    # The numbers involved are n, the coefficients (all 1), the powers of q (0 to n-1),
    # and the power of T (n-1).

    # Build the q-integer [n]_q string: e.g., (1 + q + q^2)
    q_poly_terms = ["1"]
    for i in range(1, n):
        if i == 1:
            q_poly_terms.append("q")
        else:
            q_poly_terms.append(f"q^{i}")
    
    q_poly_str = "(" + " + ".join(q_poly_terms) + ")"

    # Build the T part string: e.g., T^2
    if (n - 1) == 1:
        t_part = "T"
    else:
        t_part = f"T^{n-1}"
    
    # Print the final equation string
    print(f"nabla_q(T^{n}) = {q_poly_str} * {t_part}")

if __name__ == '__main__':
    try:
        # Prompt the user to enter an integer n
        n_input_str = input("Enter a non-negative integer n: ")
        n_input = int(n_input_str)
        display_q_derivative(n_input)
    except (ValueError, TypeError):
        print("Invalid input. Please enter a non-negative integer.")
