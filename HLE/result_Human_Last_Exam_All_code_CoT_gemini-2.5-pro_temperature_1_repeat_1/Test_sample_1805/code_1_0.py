def get_q_derivative_expression(n_str):
    """
    Computes and formats the expression for the q-derivative of T^n.
    
    Args:
        n_str (str): The integer n, provided as a string.
        
    Returns:
        str: The formatted mathematical expression.
    """
    try:
        n = int(n_str)
    except (ValueError, TypeError):
        return "Invalid input. Please provide an integer."

    # Case n = 0
    if n == 0:
        return f"∇q(T^0) = 0"

    # Determine the exponent of T
    t_exponent = n - 1

    # Determine the coefficient string for the q-integer [n]_q
    if n > 0:
        if n == 1:
            coeff_str = "1"
        else:
            # Build the polynomial string for [n]_q
            terms = ["1", "q"]
            if n > 2:
                terms.extend([f"q^{k}" for k in range(2, n)])
            coeff_str = f"({' + '.join(terms)})"
    else:  # n < 0
        m = -n
        if m == 1:
            m_terms_str = "1"
        else:
            m_terms = ["1", "q"]
            if m > 2:
                m_terms.extend([f"q^{k}" for k in range(2, m)])
            m_terms_str = f"({' + '.join(m_terms)})"
        
        # Format for negative n: [-m]_q = -q^{-m}[m]_q
        coeff_str = f"-q^{n} * {m_terms_str}"

    # Determine the string for the T part
    if t_exponent == 0:
        t_str = ""  # T^0 = 1, so it disappears
    elif t_exponent == 1:
        t_str = "T"
    else:
        t_str = f"T^{t_exponent}"

    # Assemble the final expression
    base_expr = f"∇q(T^{n})"
    if not t_str: # This happens when n=1, result is 1
        return f"{base_expr} = {coeff_str}"
    else:
        return f"{base_expr} = {coeff_str} * {t_str}"

if __name__ == '__main__':
    # Prompt the user for the integer n
    n_input = input("This script calculates ∇q(T^n).\nEnter an integer value for n: ")
    # Calculate and print the expression
    result_expression = get_q_derivative_expression(n_input)
    print(result_expression)
