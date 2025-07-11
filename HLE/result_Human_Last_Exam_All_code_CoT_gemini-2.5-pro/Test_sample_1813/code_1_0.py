import math

def compute_markov_cf(p, q):
    """
    Computes and prints the continued fraction for the generalized Markov number m_{p/q}.
    """
    print(f"Computing the continued fraction for the generalized Markov number m_{p}/{q}")
    print("-" * 60)

    # Step 1: Define the quadratic equation for the generalized Markov number.
    # The number m_p/q is the positive root of q^2*x^2 + (2p-3q)q*x - p^2 = 0.
    a = q**2
    b = (2 * p - 3 * q) * q
    c = -p**2

    print("The number is the positive root of the quadratic equation:")
    print(f"{a}x^2 + ({b})x + ({c}) = 0")
    print("-" * 60)

    # Step 2: Solve for the number using the quadratic formula.
    # The root is (-b + sqrt(b^2 - 4ac)) / 2a.
    # Let D be the discriminant.
    D = b**2 - 4 * a * c
    
    # We write the number in the form (P0 + sqrt(D)) / Q0 for the CF algorithm.
    P0 = -b
    Q0 = 2 * a
    
    print("The value of the number is:")
    print(f"m_{p}/{q} = ({P0} + sqrt({D})) / {Q0}")
    print("-" * 60)

    # Step 3: Compute the continued fraction using the standard algorithm.
    Pk = P0
    Qk = Q0
    # Use integer square root for precision
    sqrt_D_int = math.isqrt(D)

    coefficients = []
    seen_states = {}
    k = 0

    # The algorithm will terminate because the continued fraction of a quadratic
    # irrational is always periodic.
    while (Pk, Qk) not in seen_states:
        seen_states[(Pk, Qk)] = k

        # Calculate the coefficient a_k = floor((Pk + sqrt(D)) / Qk)
        ak = (Pk + sqrt_D_int) // Qk
        coefficients.append(ak)

        # Update P and Q for the next iteration
        # P_{k+1} = a_k * Q_k - P_k
        # Q_{k+1} = (D - P_{k+1}^2) / Q_k
        Pk_next = ak * Qk - Pk
        Qk_next = (D - Pk_next**2) // Qk

        Pk = Pk_next
        Qk = Qk_next
        k += 1

    # Step 4: Format and print the final continued fraction equation.
    # The number m_{4/7} is a reduced quadratic irrational, so its
    # continued fraction is purely periodic.
    start_of_period = seen_states.get((Pk, Qk), 0)
    
    pre_period = coefficients[:start_of_period]
    periodic_part = coefficients[start_of_period:]

    coeffs_str = ", ".join(map(str, periodic_part))

    print("The associated continued fraction is:")
    if not pre_period:
        # Purely periodic
        final_equation = f"m_{p}/{q} = [{coeffs_str}] (The sequence is purely periodic)"
    else:
        pre_period_str = ", ".join(map(str, pre_period))
        final_equation = f"m_{p}/{q} = [{pre_period_str}; {coeffs_str}] (The sequence starting from {periodic_part[0]} is periodic)"
    
    # Output each number in the final equation as requested.
    print(final_equation)

# Run the computation for p=4, q=7
compute_markov_cf(4, 7)