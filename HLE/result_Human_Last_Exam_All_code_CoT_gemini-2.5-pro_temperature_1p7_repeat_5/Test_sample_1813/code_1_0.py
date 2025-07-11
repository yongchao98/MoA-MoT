import math

def compute_and_display_cf(p, q):
    """
    Computes the continued fraction for p/q (for 0 < p/q < 1)
    and prints it in an expanded equation format.
    """
    if not (isinstance(p, int) and isinstance(q, int) and 0 < p < q):
        print("Error: Input must be a rational number p/q where p and q are integers and 0 < p/q < 1.")
        return

    print(f"Computing the continued fraction associated with the generalized Markov number m_{p}/{q}.")
    print(f"This is the continued fraction of the rational number {p}/{q} itself.\n")

    # Algorithm to find continued fraction coefficients
    # for p/q = 1/(k1 + 1/(k2 + ...))
    coeffs = []
    temp_p, temp_q = p, q
    while temp_p > 0:
        k = temp_q // temp_p
        coeffs.append(k)
        temp_q, temp_p = temp_p, temp_q % temp_p

    print(f"The coefficients [k1, k2, ...] for {p}/{q} are: {coeffs}\n")

    # Build and print the final equation string
    # This fulfills the requirement to "output each number in the final equation"
    equation_str = str(coeffs[-1])
    for k in reversed(coeffs[:-1]):
        equation_str = f"({k} + 1 / {equation_str})"

    print("The final equation is:")
    print(f"{p}/{q} = 1 / {equation_str}")


# Main execution for the specific case of m_{4/7}
p_val = 4
q_val = 7
compute_and_display_cf(p_val, q_val)