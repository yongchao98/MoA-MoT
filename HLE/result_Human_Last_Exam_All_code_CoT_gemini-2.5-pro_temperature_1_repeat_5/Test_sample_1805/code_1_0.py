def get_nabla_q_expression(n):
    """
    Calculates and prints the symbolic expression for the q-difference quotient of T^n.
    The formula derived is nabla_q(T^n) = [n]_q * T^(n-1), where [n]_q is the q-integer.

    Args:
        n (int): A non-negative integer representing the power of T.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: The power 'n' must be a non-negative integer.")
        return

    # Handle the case n = 0. nabla_q(T^0) = nabla_q(1) = 0.
    if n == 0:
        equation = "nabla_q(T^0) = 0"
        print(equation)
        return

    # Handle the case n = 1. nabla_q(T^1) = [1]_q * T^0 = 1 * 1 = 1.
    if n == 1:
        equation = "nabla_q(T^1) = 1"
        print(equation)
        return

    # For n > 1, we build the expression for the q-integer [n]_q.
    # [n]_q = 1 + q + q^2 + ... + q^(n-1)
    q_integer_terms = ["1"]
    if n > 1:
        q_integer_terms.append("q")
    if n > 2:
        for k in range(2, n):
            q_integer_terms.append(f"q^{k}")

    q_integer_str = f"({' + '.join(q_integer_terms)})"

    # Build the expression for the T^(n-1) part.
    if n == 2:
        t_power_str = "T"  # For T^(2-1) = T^1
    else:
        t_power_str = f"T^{n - 1}"

    # Construct the final equation string.
    # It will be of the form: nabla_q(T^n) = (1 + q + ... + q^(n-1)) * T^(n-1)
    equation = f"nabla_q(T^{n}) = {q_integer_str} * {t_power_str}"
    print(equation)


# You can change the value of n here to see the expression for different powers.
n_value = 5
get_nabla_q_expression(n_value)