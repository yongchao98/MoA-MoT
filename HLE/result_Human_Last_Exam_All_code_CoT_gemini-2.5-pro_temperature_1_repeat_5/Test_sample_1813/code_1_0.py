import math

def gcd_three(a, b, c):
    """Computes the greatest common divisor of three integers."""
    return math.gcd(math.gcd(a, b), c)

def compute_markov_cf(p, q):
    """
    Computes the continued fraction associated with the generalized Markov number m_p/q.
    
    The number is alpha = (3q - sqrt(9q^2 - 4p^2)) / (2p).
    This can be written in the form (c2 + c1*sqrt(D)) / c3.
    """
    # 1. Determine the parameters D, c1, c2, c3
    D = 9 * q**2 - 4 * p**2
    
    # Initial representation of alpha = (3q - sqrt(D)) / (2p)
    # is c1=-1, c2=3q, c3=2p
    c1 = -1
    c2 = 3 * q
    c3 = 2 * p
    
    sqrt_D_float = math.sqrt(D)
    
    # Lists to store results
    coeffs = []
    seen_states = {}

    # 2. Iteratively compute continued fraction coefficients
    while True:
        # Normalize the state (c1, c2, c3) to ensure cycle detection
        common_divisor = gcd_three(c1, c2, c3)
        c1 //= common_divisor
        c2 //= common_divisor
        c3 //= common_divisor
        
        # Make the denominator positive for a canonical representation
        if c3 < 0:
            c1, c2, c3 = -c1, -c2, -c3

        state = (c1, c2, c3)
        
        # Check if the state has been seen before
        if state in seen_states:
            start_of_period_index = seen_states[state]
            non_repeating_part = coeffs[:start_of_period_index]
            repeating_part = coeffs[start_of_period_index:]
            break
        
        seen_states[state] = len(coeffs)
        
        # Compute the integer part 'a'
        a = math.floor((c1 * sqrt_D_float + c2) / c3)
        coeffs.append(a)
        
        # Compute the next state for 1 / (xi - a)
        # xi - a = (c1*sqrt(D) + c2 - a*c3) / c3
        num_sqrt = c1
        num_const = c2 - a * c3
        
        # 1 / (xi - a) = c3 / (num_const + num_sqrt*sqrt(D))
        # Rationalize the denominator:
        # c3 * (num_const - num_sqrt*sqrt(D)) / (num_const^2 - num_sqrt^2*D)
        
        # New coefficients before normalization
        next_c1 = -c3 * num_sqrt
        next_c2 = c3 * num_const
        next_c3 = num_const**2 - num_sqrt**2 * D
        
        c1, c2, c3 = next_c1, next_c2, next_c3

    return non_repeating_part, repeating_part

# Main execution for the specific problem m_4/7
p, q = 4, 7

print(f"Computing the continued fraction associated with the generalized Markov number m_{p}/{q}.")
alpha_expr = f"({3*q} - sqrt({377})) / {2*p}"
print(f"This corresponds to the continued fraction of alpha = {alpha_expr}.")
print("-" * 20)

non_repeating, repeating = compute_markov_cf(p, q)

# Format the output string for the continued fraction
non_repeating_str = ', '.join(map(str, non_repeating[1:]))
repeating_str = ', '.join(map(str, repeating))
cf_string = f"[{non_repeating[0]}; "
if non_repeating_str:
    cf_string += f"{non_repeating_str}, "
cf_string += f"overline({repeating_str})]"

print(f"The resulting continued fraction is: {cf_string}\n")

print("The numbers (coefficients) in the final continued fraction are:")
print(f"Non-repeating part: {non_repeating}")
print(f"Repeating part: {repeating}")

print("\nDetailed coefficients:")
for i, val in enumerate(non_repeating):
    print(f"a_{i} = {val}")

print("The repeating block starts after the non-repeating part:")
for i, val in enumerate(repeating):
    print(f"a_{i + len(non_repeating)} = {val} (repeats)")
