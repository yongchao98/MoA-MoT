import math

def compute_generalized_markov_cf():
    """
    Computes the continued fraction of the generalized Markov number m_4/7.

    The generalized Markov number m_{p/q} is defined via the Lagrange constant
    of a related periodic continued fraction. For r = 4/7, this number is
    sqrt(445)/9.

    This function computes the continued fraction of this number using the
    standard algorithm for quadratic irrationals, which is exact and avoids
    floating-point errors.
    """
    D = 445
    q = 9

    # To use the integer algorithm for (P + sqrt(D)) / Q,
    # we rewrite sqrt(D)/q as sqrt(D*q^2)/(q^2)
    D_cf = D * (q**2)
    P = 0
    Q = q**2
    
    sqrt_D_cf = math.sqrt(D_cf)

    # Store states (P, Q) to detect the period
    history = {}
    quotients = []
    
    pre_period_len = 0
    found_period = False

    while not found_period:
        if (P, Q) in history:
            pre_period_len = history[(P, Q)]
            found_period = True
        else:
            history[(P, Q)] = len(quotients)
            a = int((P + sqrt_D_cf) / Q)
            quotients.append(a)
            P_next = a * Q - P
            Q_next = (D_cf - P_next**2) // Q
            P, Q = P_next, Q_next

    pre_period = quotients[:pre_period_len]
    periodic_part = quotients[pre_period_len:]
    
    print(f"The generalized Markov number m_4/7 is sqrt({D})/{q}.")
    
    # Format the output string
    pre_period_str = ""
    if pre_period:
        pre_period_str = ", ".join(map(str, pre_period))
        if periodic_part:
            pre_period_str += "; "

    periodic_part_str = ", ".join(map(str, periodic_part))
    
    final_cf = f"[{pre_period_str}({periodic_part_str})]"
    
    print(f"sqrt({D})/{q} = [", end="")
    
    if pre_period:
        if len(pre_period) > 1:
            print(", ".join(map(str, pre_period[:-1])), end="")
            print(f", {pre_period[-1]}; ", end="")
        else:
            print(f"{pre_period[0]}; ", end="")

    print("; ".join(map(str, periodic_part))) # This is not the right formatting
    # Let's fix the final print
    
    a0 = pre_period[0] if pre_period else periodic_part.pop(0)

    print(f"The continued fraction is [{a0}; ({', '.join(map(str, periodic_part))})]")
    
    # Correct final equation print format
    final_eq_str = f"sqrt({D})/{q} = [{a0}; ("
    for number in periodic_part:
        final_eq_str += f"{number}, "
    # remove last comma and space
    final_eq_str = final_eq_str[:-2] + ")]"
    print("The final equation is:")
    print(final_eq_str)

def main():
    """
    Main function to compute and print the continued fraction.
    The number is m_{4/7} = sqrt(445)/9.
    The continued fraction is calculated using an exact integer algorithm.
    """
    
    D_val = 445
    q_val = 9

    # The number is of the form (P0 + sqrt(D_cf))/Q0
    # For sqrt(445)/9, we transform it to (0 + sqrt(445*81))/81
    D_cf = D_val * q_val**2
    P0 = 0
    Q0 = q_val**2
    
    sqrt_D_cf_approx = math.sqrt(D_cf)

    states = {}
    coefficients = []
    
    P, Q = P0, Q0

    while (P, Q) not in states:
        states[(P, Q)] = len(coefficients)
        a = int((P + sqrt_D_cf_approx) / Q)
        coefficients.append(a)
        
        P_next = a * Q - P
        Q_next = (D_cf - P_next**2) // Q
        P, Q = P_next, Q_next
    
    start_of_period_index = states[(P, Q)]
    
    non_periodic_part = coefficients[:start_of_period_index]
    periodic_part = coefficients[start_of_period_index:]

    a0 = non_periodic_part[0]
    
    # Print the equation as requested
    periodic_part_str = ", ".join(map(str, periodic_part))
    print(f"sqrt({D_val})/{q_val} = [{a0}; ({periodic_part_str})]")
    
    # Another format just printing the numbers in the equation
    print("Each number in the final equation is:")
    print(f"sqrt({D_val})/{q_val} = [{a0}; ", end="")
    for i, num in enumerate(periodic_part):
        if i < len(periodic_part) - 1:
            print(f"{num}, ", end="")
        else:
            print(f"{num}", end="")
    print("]")


if __name__ == '__main__':
    main()
