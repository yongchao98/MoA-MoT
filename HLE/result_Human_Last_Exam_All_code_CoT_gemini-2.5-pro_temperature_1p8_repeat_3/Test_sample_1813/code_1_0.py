import math

def generalized_markov_cf(p_frac, q_frac):
    """
    Computes the continued fraction of the generalized Markov number m_{p/q}.
    """
    # Step 1: Find modular inverse p' such that p*p' === 1 (mod q)
    try:
        # This function is available in Python 3.8+
        pp_inv = pow(p_frac, -1, q_frac)
    except AttributeError:
        # Fallback for older Python versions
        def egcd(a, b):
            if a == 0:
                return (b, 0, 1)
            g, y, x = egcd(b % a, a)
            return (g, x - (b // a) * y, y)

        def modinv(a, m):
            g, x, y = egcd(a, m)
            if g != 1:
                raise Exception('modular inverse does not exist')
            return x % m
        pp_inv = modinv(p_frac, q_frac)
    
    # Step 2: Use Frobenius's formula to get the quadratic equation Az^2+Bz+C=0
    # The coefficients are simplified after multiplying by q_frac.
    A = q_frac**2
    B = -(3 * q_frac**2 - 2 * pp_inv * q_frac)
    C = (3 * pp_inv * q_frac - p_frac * q_frac - (pp_inv**2 + 1))
    
    # The number is the larger root: z = (-B + sqrt(B^2-4AC))/(2A)
    # Match this to the PQa algorithm's input: (P0 + sqrt(D))/Q0
    D_sqrt = B**2 - 4*A*C
    P0 = -B
    Q0 = 2*A

    # Step 3: PQa Algorithm to compute the continued fraction
    coeffs = []
    seen_states = {}
    
    pi = P0
    qi = Q0
    
    # Loop until a state (pi, qi) repeats
    while (pi, qi) not in seen_states:
        seen_states[(pi, qi)] = len(coeffs)
        
        # Calculate the next coefficient
        ai = int((pi + math.sqrt(D_sqrt)) / qi)
        coeffs.append(ai)
        
        # Calculate the next state
        pi_next = ai * qi - pi
        qi_next = (D_sqrt - pi_next**2) // qi
        
        pi = pi_next
        qi = qi_next
        
    # The period starts where the state first repeated
    period_start_index = seen_states[(pi, qi)]
    non_periodic_part = coeffs[:period_start_index]
    periodic_part = coeffs[period_start_index:]
    
    return non_periodic_part, periodic_part

# Main execution
p, q = 4, 7
non_periodic, periodic = generalized_markov_cf(p, q)

# Print the final results as requested
print(f"The continued fraction associated with the generalized Markov number m_{p}/{q} is:")

# Constructing the output string
result_string = f"[{non_periodic[0]}; ("
result_string += ", ".join(map(str, periodic))
result_string += ")]"
print(result_string)

print("\nThe coefficients of the continued fraction are:")
full_coeffs = non_periodic + periodic
print(f"a_0 = {non_periodic[0]}")
# To fulfill the requirement of printing each coefficient of the final equation
for i in range(len(periodic)):
    # Index relative to start of CF
    idx = i + len(non_periodic)
    print(f"a_{idx} = {periodic[i]}")
