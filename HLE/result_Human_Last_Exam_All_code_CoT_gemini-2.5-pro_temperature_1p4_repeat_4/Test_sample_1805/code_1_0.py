def calculate_and_print_q_nabla(n: int):
    """
    Calculates and prints the symbolic expression for nabla_q(T^n) for a non-negative integer n.

    The q-derivative nabla_q(T^n) is given by [n]_q * T^(n-1), where [n]_q
    is the q-integer. For non-negative n, [n]_q can be expanded as the
    polynomial 1 + q + q^2 + ... + q^(n-1).
    """
    if not isinstance(n, int) or n < 0:
        print("This script is designed for non-negative integers n.")
        return

    # Handle the case n=0. nabla_q(T^0) = nabla_q(1) = 0.
    if n == 0:
        print(f"nabla_q(T^{n}) = 0")
        return

    # Handle the case n=1. nabla_q(T^1) = [1]_q * T^0 = 1 * 1 = 1.
    if n == 1:
        print(f"nabla_q(T^{n}) = 1")
        return

    # For n > 1, build the string for the q-polynomial [n]_q
    q_poly_parts = []
    for k in range(n):
        if k == 0:
            q_poly_parts.append("1")
        elif k == 1:
            q_poly_parts.append("q")
        else:
            # Append q^k for k >= 2
            q_poly_parts.append(f"q^{k}")
    
    # Join the parts with ' + ' and enclose in parentheses
    q_poly_str = f"({' + '.join(q_poly_parts)})"

    # Build the string for the T part, T^(n-1)
    power = n - 1
    if power == 1:
        t_str = " * T"
    else:
        t_str = f" * T^{power}"

    # Print the final formatted equation
    print(f"nabla_q(T^{n}) = {q_poly_str}{t_str}")


if __name__ == '__main__':
    # You can change the value of n here to see the result for different integers.
    n_example = 5
    calculate_and_print_q_nabla(n_example)
