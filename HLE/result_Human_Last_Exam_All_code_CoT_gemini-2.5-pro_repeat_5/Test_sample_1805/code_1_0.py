def generate_nabla_q_expression():
    """
    Asks the user for a non-negative integer n and prints the symbolic
    expression for the q-derivative of T^n. The context Z[T] implies
    n is a non-negative integer for polynomials.
    """
    try:
        # Prompt the user for input
        n_input = input("Please provide a non-negative integer n for T^n: ")
        n = int(n_input)
    except (ValueError, TypeError):
        print("Error: Invalid input. Please enter an integer.")
        return

    # In the context of polynomials Z[T], n is typically non-negative.
    if n < 0:
        print("Error: This calculation is for non-negative integers as exponents in polynomials.")
        return

    # Handle the special case n = 0
    if n == 0:
        # The q-derivative of a constant T^0 = 1 is 0.
        # This aligns with the formula [0]_q = (q^0 - 1)/(q - 1) = 0.
        print("The final equation is: nabla_q(T^0) = 0")
        return

    # Build the string for the q-analogue [n]_q = 1 + q + ... + q^(n-1)
    q_terms = []
    # Add the first term "1" (for q^0)
    q_terms.append("1")
    if n > 1:
        # Add the second term "q" (for q^1)
        q_terms.append("q")
    if n > 2:
        # Add subsequent terms "q^i"
        for i in range(2, n):
            q_terms.append(f"q^{i}")

    if len(q_terms) > 1:
        # Enclose in parentheses if it's a sum of multiple terms
        q_part_str = "(" + " + ".join(q_terms) + ")"
    else:
        # Just "1" if n=1
        q_part_str = q_terms[0]

    # Build the string for the T part, T^(n-1)
    t_power = n - 1
    t_part_str = ""
    if t_power == 1:
        t_part_str = "T"
    elif t_power > 1:
        t_part_str = f"T^{t_power}"
    # If t_power is 0, t_part_str remains "", representing T^0 = 1.

    # Assemble the final equation string: nabla_q(T^n) = ...
    equation_lhs = f"nabla_q(T^{n})"

    # Combine the q-analogue part and the T part
    if t_part_str:  # This is true if n > 1
        equation_rhs = f"{q_part_str} * {t_part_str}"
    else:  # This is true if n = 1 (t_power is 0)
        equation_rhs = q_part_str

    final_equation = f"{equation_lhs} = {equation_rhs}"

    print("The final equation is:")
    print(final_equation)

if __name__ == '__main__':
    generate_nabla_q_expression()