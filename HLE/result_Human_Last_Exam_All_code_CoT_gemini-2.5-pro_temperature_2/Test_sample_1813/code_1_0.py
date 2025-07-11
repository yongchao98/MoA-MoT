import math

def solve_task():
    """
    Computes the continued fraction associated with the generalized Markov number m_{4/7}.
    """
    # Step 1 & 2: Set up the equation for the Markov number m_{4/7}
    p, q = 4, 7
    # We solve 4*q_prime - 7*p_prime = 1. A solution is p_prime=1, q_prime=2.
    p_prime, q_prime = 1, 2
    
    # Quadratic equation is x^2 - 3*p*q_prime*x + (p^2 + q_prime^2) = 0
    # x^2 - 3*4*2*x + (4^2 + 2^2) = 0
    # x^2 - 24x + 20 = 0
    
    # Step 3: Solve the quadratic equation to find m_{4/7}
    a, b, c = 1, -24, 20
    discriminant = b**2 - 4*a*c  # 496
    
    # The number m_{4/7} is the larger root: (24 + sqrt(496))/2 = 12 + sqrt(124)
    # sqrt(496) = sqrt(16 * 31) = 4 * sqrt(31)
    # The number is 12 + 2*sqrt(31)
    
    # Step 4: Compute the continued fraction of 12 + sqrt(124)
    # We use the standard algorithm for (A + sqrt(D)) / B
    
    # Initial state for xi_0 = (12 + sqrt(124)) / 1
    D = 124
    A = 12
    B = 1
    
    coeffs = []
    seen_states = {}

    sqrt_D_int = math.isqrt(D)

    k = 0
    while (A, B) not in seen_states:
        seen_states[(A, B)] = k
        
        # Calculate the integer part of the current number
        ak = (A + sqrt_D_int) // B
        coeffs.append(ak)
        
        # Update A and B for the next iteration
        next_A = ak * B - A
        next_B = (D - next_A**2) // B
        
        A = next_A
        B = next_B
        k += 1

    pre_period_len = seen_states[(A, B)]
    a0 = coeffs[0]
    pre_period = coeffs[1:pre_period_len] # Should be empty in this case
    period = coeffs[pre_period_len:]

    print(f"The generalized Markov number m_{p}/{q} is the larger root of the equation x^2 - {3*p*q_prime}x + {p**2 + q_prime**2} = 0.")
    print(f"The value is m_{p}/{q} = 12 + 2*sqrt(31).")
    
    period_str = ", ".join(map(str, period))
    if pre_period:
        pre_period_str = ", ".join(map(str, pre_period))
        print(f"Its continued fraction is: [{a0}; {pre_period_str}, ({period_str})]")
    else:
        print(f"Its continued fraction is: [{a0}; ({period_str})]")

    print("\nThe numbers in the final continued fraction equation are:")
    print(f"a_0 = {coeffs[0]}")
    for i in range(1, len(coeffs)):
        print(f"a_{i} = {coeffs[i]}")

solve_task()
<<<[23; (7, 2, 1, 1, 1, 3, 1, 4, 1, 3, 1, 1, 1, 2, 7, 22)]>>>