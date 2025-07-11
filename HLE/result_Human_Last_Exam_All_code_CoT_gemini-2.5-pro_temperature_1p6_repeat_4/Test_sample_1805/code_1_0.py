def generate_nabla_q_expression(n):
    """
    Generates and prints the mathematical expression for the q-derivative of T^n.

    Args:
        n (int): The integer exponent, must be non-negative.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: The exponent n must be a non-negative integer.")
        return

    # Case n = 0: The derivative of a constant T^0 = 1 is 0.
    if n == 0:
        print("nabla_q(T^0) = 0")
        return

    # Case n = 1: nabla_q(T^1) = [1]_q * T^0 = 1
    if n == 1:
        print("nabla_q(T^1) = 1")
        return

    # General case for n > 1
    # Build the q-integer [n]_q part: (1 + q + q^2 + ... + q^(n-1))
    q_terms = []
    # i = 0 gives q^0 = 1
    q_terms.append("1")
    # i = 1 gives q^1 = q
    if n > 1:
        q_terms.append("q")
    # i > 1 gives q^i
    for i in range(2, n):
        q_terms.append(f"q^{i}")
    
    coeff_str = f"({ ' + '.join(q_terms) })"

    # Build the T part: T^(n-1)
    t_exp = n - 1
    if t_exp == 1:
        t_str = "T"
    else:
        t_str = f"T^{t_exp}"
    
    # Combine and print the final expression
    print(f"nabla_q(T^{n}) = {coeff_str}{t_str}")


# You can change this value to see the expression for different n.
n_value = 5
generate_nabla_q_expression(n_value)
