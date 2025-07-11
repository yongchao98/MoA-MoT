def calculate_q_derivative_of_T_n(n):
    """
    Calculates and prints the q-derivative of T^n.
    The result is expressed as nabla_q(T^n) = [n]_q * T^(n-1),
    where [n]_q is the q-integer 1 + q + q^2 + ... + q^(n-1).
    """
    if not isinstance(n, int) or n < 0:
        print("Error: n must be a non-negative integer.")
        return

    # Handle the case n=0
    if n == 0:
        # nabla_q(T^0) = nabla_q(1) = (1 - 1) / ((q-1)T) = 0
        print("∇_q(T^0) = 0")
        return

    # Handle the case n=1
    if n == 1:
        # nabla_q(T^1) = (qT - T) / (qT - T) = 1
        print("∇_q(T^1) = 1")
        return

    # General case for n > 1
    
    # Build the q-integer polynomial string: [n]_q
    q_poly_terms = []
    # Term for i=0 is q^0 = 1
    q_poly_terms.append("1")
    # Term for i=1 is q^1 = q
    q_poly_terms.append("q")
    # Terms for i > 1 are q^i
    for i in range(2, n):
        q_poly_terms.append(f"q^{i}")
    
    q_poly_str = " + ".join(q_poly_terms)

    # Build the T part string: T^(n-1)
    t_power = n - 1
    if t_power == 1:
        t_part = "T"
    else:
        t_part = f"T^{t_power}"

    # Construct the final equation string
    # The result shows each individual term in the q-integer expansion
    final_equation = f"∇_q(T^{n}) = ({q_poly_str}) * {t_part}"
    
    print(final_equation)

# --- Main execution ---
# Let's calculate the result for a specific integer, e.g., n=5
n_value = 5
calculate_q_derivative_of_T_n(n_value)