def generate_nabla_q_formula():
    """
    Prompts the user for an integer n and prints the formula for nabla_q(T^n).
    """
    try:
        n_str = input("Enter a non-negative integer n: ")
        n = int(n_str)
        if n < 0:
            print("Error: n must be a non-negative integer.")
            return
    except ValueError:
        print(f"Error: Invalid input '{n_str}'. Please enter an integer.")
        return

    # The left-hand side of the equation
    lhs = f"nabla_q(T^{n})"

    # Handle the case n=0, where nabla_q(constant) = 0
    if n == 0:
        print(f"The equation for n={n} is:")
        print(f"nabla_q(T^0) = 0")
        return

    # Build the q-analogue sum part: [n]_q = (1 + q + ... + q^(n-1))
    if n == 1:
        q_sum_str = "1"
    else:
        terms = ["1"]
        # Loop to create terms q, q^2, ..., q^(n-1)
        for i in range(1, n):
            if i == 1:
                terms.append("q")
            else:
                terms.append(f"q^{i}")
        q_sum_str = f"({' + '.join(terms)})"

    # Build the T term part: T^(n-1)
    if n == 1:
        # T^0 which is 1
        t_term_str = "T^0"
    elif n == 2:
        # T^1 which is T
        t_term_str = "T"
    else:
        t_term_str = f"T^{n-1}"
    
    # Assemble the right-hand side of the equation
    rhs = f"{q_sum_str} * {t_term_str} * nabla_q(T)"

    print(f"\nThe formula for n={n} is:")
    print(f"{lhs} = {rhs}")
    
    # Add a simplified version for n=1 for clarity
    if n == 1:
        print("Which simplifies to:")
        print("nabla_q(T) = nabla_q(T)")


if __name__ == '__main__':
    generate_nabla_q_formula()