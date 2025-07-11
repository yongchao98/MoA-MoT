def generate_q_nabla_expression(n: int):
    """
    Calculates and prints the symbolic expression for the q-difference quotient of T^n.

    The q-difference quotient (or Jackson derivative) nabla_q of a monomial T^n is given by:
    nabla_q(T^n) = [n]_q * T^(n-1)
    where [n]_q is the q-integer, defined as (1 + q + q^2 + ... + q^(n-1)).

    Args:
        n: A non-negative integer for the power of T.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: The input 'n' must be a non-negative integer.")
        return

    # Handle the n = 0 case: nabla_q(T^0) = nabla_q(1) = 0
    if n == 0:
        # The equation is nabla_q(T^0) = 0
        print(f"nabla_q(T^{0}) = 0")
        return

    # Build the q-integer polynomial part: [n]_q = 1 + q + ... + q^(n-1)
    q_terms = []
    # The sum goes from i=0 to n-1
    for i in range(n):
        if i == 0:
            q_terms.append("1")
        elif i == 1:
            q_terms.append("q")
        else:
            q_terms.append(f"q^{i}")
    q_polynomial_str = " + ".join(q_terms)

    # Build the T term: T^(n-1)
    t_power = n - 1
    t_term_str = ""
    if t_power == 1:
        # For T^1, just show T
        t_term_str = " * T"
    elif t_power > 1:
        # For T^k where k > 1
        t_term_str = f" * T^{t_power}"
    # If t_power is 0, t_term_str remains empty, as T^0 = 1

    # Combine the parts into the final equation string
    # For n=1, [1]_q = 1 and T^0=1, so the result is 1. Parentheses are not needed.
    if n == 1:
        result_equation = f"nabla_q(T^{n}) = {q_polynomial_str}"
    # For n > 1, use parentheses for clarity
    else:
        result_equation = f"nabla_q(T^{n}) = ({q_polynomial_str}){t_term_str}"

    print(result_equation)


# --- Main execution ---
# Set the integer 'n' for which to calculate nabla_q(T^n).
# You can change this value to any non-negative integer.
n_value = 5

# Generate and print the expression for the chosen n.
generate_q_nabla_expression(n_value)