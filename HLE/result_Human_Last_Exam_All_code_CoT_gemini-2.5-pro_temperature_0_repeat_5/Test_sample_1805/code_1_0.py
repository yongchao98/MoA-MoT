def q_nabla_expression(n):
    """
    This function generates and prints the mathematical expression for nabla_q(T^n).

    Args:
        n (int): A non-negative integer.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: n must be a non-negative integer.")
        return

    # Handle the case n=0
    if n == 0:
        # nabla_q(T^0) = nabla_q(1) = (1-1)/((q-1)T) = 0
        print("nabla_q(T^0) = 0")
        return

    # Build the q-integer part: [n]_q = 1 + q + q^2 + ... + q^(n-1)
    q_terms = []
    for i in range(n):
        if i == 0:
            q_terms.append("1")
        elif i == 1:
            q_terms.append("q")
        else:
            q_terms.append(f"q^{i}")
    
    # For n=1, [1]_q = 1. For n>1, it's a sum.
    if n == 1:
        q_poly_str = "1"
    else:
        q_poly_str = f"({' + '.join(q_terms)})"

    # Build the T part: T^(n-1)
    if n - 1 == 0:
        # For n=1, T^0 = 1, so no T term is printed
        t_term_str = ""
    elif n - 1 == 1:
        t_term_str = " * T"
    else:
        t_term_str = f" * T^{n-1}"

    # Construct the final equation string
    final_equation = f"nabla_q(T^{n}) = {q_poly_str}{t_term_str}"

    print(final_equation)

# --- User-editable section ---
# You can change the value of n to see the result for different integers.
integer_n = 5
# --- End of user-editable section ---

print(f"Calculating the expression for n = {integer_n}:")
q_nabla_expression(integer_n)