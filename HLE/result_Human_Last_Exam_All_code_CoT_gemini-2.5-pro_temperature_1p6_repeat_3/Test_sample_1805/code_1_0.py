def q_derivative_of_power(n: int):
    """
    Calculates and prints the formula for the q-derivative of T^n.

    The q-derivative is given by the formula:
    ∇_q(T^n) = Σ_{k=0}^{n-1} q^k * T^{n-1}

    Args:
        n: A non-negative integer.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: n must be a non-negative integer.")
        return

    # Case n=0: ∇_q(T^0) = ∇_q(1) = 0
    if n == 0:
        print(f"∇_q(T^{n}) = 0")
        return

    # Case n=1: ∇_q(T^1) = T^0 = 1
    if n == 1:
        print(f"∇_q(T^{n}) = 1")
        return

    # General case for n > 1
    t_exponent = n - 1
    
    # Build the string for the T part of the expression
    if t_exponent == 1:
        t_str = "T"
    else:
        t_str = f"T^{t_exponent}"

    # Build the list of terms in the sum
    terms = []
    # First term (k=0), which is q^0 * T^(n-1) = T^(n-1)
    terms.append(t_str)

    # Subsequent terms (k=1 to n-1)
    for k in range(1, n):
        if k == 1:
            q_str = "q"
        else:
            q_str = f"q^{k}"
        
        # We output each number in the equation, as requested
        terms.append(f"{q_str} * {t_str}")

    # Join all the terms with " + " to form the right-hand side of the equation
    rhs = " + ".join(terms)
    
    print(f"∇_q(T^{n}) = {rhs}")


# Example: Calculate ∇_q(T^n) for a specific integer n.
# Let's demonstrate with n = 5.
n_value = 5
q_derivative_of_power(n_value)