import math

def get_cf_rational(num, den):
    """Computes the simple continued fraction coefficients for a rational number."""
    coeffs = []
    if den == 0:
        raise ValueError("Denominator cannot be zero.")
    while den != 0:
        a = num // den
        coeffs.append(a)
        num, den = den, num % den
    return coeffs

def get_last_two_convergents(coeffs):
    """
    Computes the last two convergents for a sequence of CF coefficients [c_1; c_2, ...].
    Returns ((p_k, q_k), (p_{k-1}, q_{k-1})).
    """
    if not coeffs:
        # Represents an empty CF, which is conventionally 0.
        # Convergents are p0/q0=0/1, p-1/q-1=1/0
        return (0, 1), (1, 0)
    
    # Initialization for p_k/q_k = [c_1; c_2, ..., c_k]
    p_km1, q_km1 = 1, 0  # Zeroth convergent p_0/q_0
    p_k, q_k = coeffs[0], 1  # First convergent p_1/q_1
    
    if len(coeffs) == 1:
        return (p_k, q_k), (p_km1, q_km1)
        
    for i in range(1, len(coeffs)):
        a = coeffs[i]
        p_next = a * p_k + p_km1
        q_next = a * q_k + q_km1
        
        p_km1, q_km1 = p_k, q_k
        p_k, q_k = p_next, q_next
    
    return (p_k, q_k), (p_km1, q_km1)

def simplify_sqrt(n):
    """Simplifies sqrt(n) to a * sqrt(b), returning (a, b)."""
    if n < 0:
        raise ValueError("Cannot take the square root of a negative number.")
    i = int(math.sqrt(n))
    while i > 1:
        if n % (i * i) == 0:
            return i, n // (i * i)
        i -= 1
    return 1, n

def main():
    """
    Main function to compute the continued fraction associated with m_{4/7}.
    """
    p, q = 4, 7

    # Step 1: Get CF of the rational number
    cf_rat = get_cf_rational(p, q)
    a0 = cf_rat[0]
    seq = cf_rat[1:]
    
    print(f"The continued fraction of {p}/{q} is [{cf_rat[0]}; {', '.join(map(str, cf_rat[1:]))}].")

    # Step 2: Construct the periodic part for xi
    period = seq + seq[::-1]
    period_str = ", ".join(map(str, period))
    print(f"The associated periodic part is ({period_str}).")
    print(f"The full continued fraction is [{a0}; {period_str}, ...].")

    # Step 3: Compute the value of the periodic continued fraction
    (p_k, q_k), (p_km1, q_km1) = get_last_two_convergents(period)
    
    # The value y = [overline(period)] satisfies the equation:
    # y = (p_k * y + p_km1) / (q_k * y + q_km1)
    # This leads to the quadratic equation:
    # q_k * y^2 + (q_km1 - p_k) * y - p_km1 = 0
    A = q_k
    B = q_km1 - p_k
    C = -p_km1
    
    # Solve the quadratic equation for y using the quadratic formula.
    # We take the positive root since the coefficients are positive.
    # y = (-B + sqrt(B^2 - 4AC)) / 2A
    discriminant = B**2 - 4*A*C
    
    # Simplify the radical part of the solution
    s_coeff, s_root = simplify_sqrt(discriminant)

    # Simplify the expression for y
    y_num_const = -B
    y_num_sqrt_coeff = s_coeff
    y_den = 2*A
    common_divisor_y = math.gcd(math.gcd(y_num_const, y_num_sqrt_coeff), y_den)
    y_const = y_num_const // common_divisor_y
    y_sqrt_coeff = y_num_sqrt_coeff // common_divisor_y
    y_den = y_den // common_divisor_y
    
    # The number we seek is xi = a0 + 1/y. Since a0=0, xi = 1/y.
    # xi = y_den / (y_const + y_sqrt_coeff * sqrt(s_root))
    
    # Rationalize the denominator to get xi in a standard form
    xi_num_const = y_den * y_const
    xi_num_sqrt_coeff = -y_den * y_sqrt_coeff
    xi_den = y_const**2 - (y_sqrt_coeff**2 * s_root)
    
    # Final simplification of xi
    common_divisor_xi = math.gcd(math.gcd(xi_num_const, xi_num_sqrt_coeff), xi_den)
    final_const = xi_num_const // common_divisor_xi
    final_sqrt_coeff = xi_num_sqrt_coeff // common_divisor_xi
    final_den = xi_den // common_divisor_xi
    
    # Ensure the denominator is positive for standard form
    if final_den < 0:
        final_const = -final_const
        final_sqrt_coeff = -final_sqrt_coeff
        final_den = -final_den

    # Format the output expression for readability
    if final_sqrt_coeff == 1:
        sqrt_part = f"sqrt({s_root})"
    elif final_sqrt_coeff == -1:
        sqrt_part = f"-sqrt({s_root})"
    else:
        sqrt_part = f"{final_sqrt_coeff}*sqrt({s_root})"

    if final_const < 0:
        # Put positive sqrt term first
        expr_str = f"({sqrt_part} - {abs(final_const)}) / {final_den}"
    elif final_const > 0:
        expr_str = f"({final_const} + {sqrt_part}) / {final_den}"
    else: # final_const is 0
        expr_str = f"{sqrt_part} / {final_den}"

    print(f"\nThis continued fraction evaluates to the number: {expr_str}")

    # Output each number in the final equation as requested
    print("\nThe numbers appearing in the final equation are:")
    all_numbers = [a0] + period + [s_root, abs(final_const), final_den]
    for num in all_numbers:
      print(num)

if __name__ == '__main__':
    main()