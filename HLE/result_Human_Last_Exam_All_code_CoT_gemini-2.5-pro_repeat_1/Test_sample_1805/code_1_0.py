def calculate_q_derivative_of_monomial(n: int):
    """
    Calculates and prints the q-derivative of T^n for a non-negative integer n.

    The q-derivative of T^n is [n]_q * T^(n-1), where [n]_q is the q-integer.
    [n]_q = 1 + q + q^2 + ... + q^(n-1).
    """
    # The context Z[T] implies T is a variable in a polynomial,
    # so n must be a non-negative integer.
    if not isinstance(n, int) or n < 0:
        print("Error: This function is defined for non-negative integers n.")
        return

    # Left-hand side of the equation
    lhs = f"nabla_q(T^{n})"

    # Handle the special case n = 0
    # nabla_q(T^0) = nabla_q(1) = 0
    if n == 0:
        rhs = "0"
        print(f"For n = {n}:")
        print(f"{lhs} = {rhs}")
        return
        
    # Handle the special case n = 1
    # nabla_q(T^1) = [1]_q * T^0 = 1 * 1 = 1
    if n == 1:
        rhs = "1"
        print(f"For n = {n}:")
        print(f"{lhs} = {rhs}")
        return

    # General case for n > 1
    
    # Build the string for the q-integer [n]_q
    q_integer_terms = []
    # Term for k=0
    q_integer_terms.append("1")
    # Term for k=1
    q_integer_terms.append("q")
    # Terms for k > 1
    for k in range(2, n):
        q_integer_terms.append(f"q^{k}")
    
    q_integer_str = " + ".join(q_integer_terms)

    # Build the string for the T part
    power_of_T = n - 1
    if power_of_T == 1:
        t_str = "T"
    else:
        t_str = f"T^{power_of_T}"

    # Combine to form the right-hand side of the equation
    rhs = f"({q_integer_str}) * {t_str}"

    print(f"For n = {n}:")
    print(f"{lhs} = {rhs}")

# --- Main execution ---
# You can change this value to see the result for different integers.
n_value = 5
calculate_q_derivative_of_monomial(n_value)
