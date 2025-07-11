def calculate_q_derivative_expression(n_str):
    """
    Calculates and prints the expression for the q-derivative of T^n.

    The q-derivative, or Jackson derivative, of T^n is given by:
    ∇_q(T^n) = [n]_q * T^(n-1)
    where [n]_q is the q-analog of n, equal to (1 + q + q^2 + ... + q^(n-1)).

    This function constructs a string representation of this formula for a given integer n.
    """
    try:
        n = int(n_str)
    except ValueError:
        print("Error: Please enter a valid integer for n.")
        return

    # Handle n = 0
    if n == 0:
        # ∇_q(T^0) = ∇_q(1) = (1 - 1) / ((q-1)T) = 0
        print("∇_q(T^0) = 0")
        return

    # Handle n = 1
    if n == 1:
        # [1]_q = 1, T^(1-1) = T^0 = 1
        print("∇_q(T^1) = T^0 = 1")
        return

    # Build the q-analog part, [n]_q
    q_terms = []
    if n > 0:
        q_terms.append("1")
        if n > 1:
            q_terms.append("q")
        for i in range(2, n):
            q_terms.append(f"q^{i}")
    # For negative n, the formula is more complex, but we can write the formal sum
    elif n < 0:
         # [n]_q = -q^n * [-n]_(1/q)
         # We'll just present the formal result for simplicity here.
         print(f"For n < 0, the expression is more complex. Formally, ∇_q(T^{n}) = [n]_q * T^(n-1)")
         return


    q_analog_str = " + ".join(q_terms)

    # Build the T part, T^(n-1)
    if (n - 1) == 0:
        t_power_str = "1" # T^0 = 1
    elif (n - 1) == 1:
        t_power_str = "T" # T^1 = T
    else:
        t_power_str = f"T^{n-1}"

    # Combine into the final expression
    if t_power_str == "1":
         final_expr = f"∇_q(T^{n}) = {q_analog_str}"
    else:
         final_expr = f"∇_q(T^{n}) = ({q_analog_str}) * {t_power_str}"

    print(final_expr)


if __name__ == '__main__':
    # Get user input for n
    n_input = input("Enter an integer value for n: ")
    calculate_q_derivative_expression(n_input)
