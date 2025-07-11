import math

def gcd_three(a, b, c):
    """Computes the greatest common divisor of three integers."""
    return math.gcd(math.gcd(a, b), c)

def get_continued_fraction_for_quadratic_irrational(p, q):
    """
    Computes the continued fraction for the quadratic irrational associated with m_{p/q}.
    The number is alpha = (3 - sqrt(9 - 4*(p/q)**2)) / (2*(p/q)).
    This simplifies to (3*q - sqrt(9*q**2 - 4*p**2)) / (2*p).
    We represent a number (A + B*sqrt(D)) / C as a tuple (A, B, C).
    """
    # D = 9*q^2 - 4*p^2
    D = 9 * q**2 - 4 * p**2
    
    # Initial number: (3q - sqrt(D)) / (2p)
    # (A, B, C) represents (A + B*sqrt(D)) / C
    # So, A = 3q, B = -1, C = 2p
    A = 3 * q
    B = -1
    C = 2 * p
    
    # Simplify the initial representation
    common_divisor = gcd_three(A, B, C)
    A //= common_divisor
    B //= common_divisor
    C //= common_divisor

    # Store seen states (A, B, C) to detect cycles
    seen_states = {}
    coeffs = []
    
    # Loop to find coefficients and detect the cycle
    while (A, B, C) not in seen_states:
        seen_states[(A, B, C)] = len(coeffs)
        
        # Calculate the integer part a_i
        # a_i = floor((A + B*sqrt(D)) / C)
        # Note: B can be negative
        a_i = math.floor((A + B * math.sqrt(D)) / C)
        coeffs.append(a_i)
        
        # Calculate the next state for alpha_{i+1} = 1 / (alpha_i - a_i)
        # alpha_i - a_i = (A - a_i*C + B*sqrt(D)) / C
        # Let A' = A - a_i*C
        A_prime = A - a_i * C
        
        # 1 / (alpha_i - a_i) = C / (A' + B*sqrt(D))
        # Rationalize the denominator: C * (A' - B*sqrt(D)) / (A'^2 - B^2*D)
        next_A_num = C * A_prime
        next_B_num = -C * B
        next_C_den = A_prime**2 - B**2 * D
        
        # Update A, B, C for the next iteration
        A = next_A_num
        B = next_B_num
        C = next_C_den
        
        # Simplify the new fraction by dividing by GCD
        common_divisor = gcd_three(A, B, C)
        A //= common_divisor
        B //= common_divisor
        C //= common_divisor

    # Cycle detected
    start_of_cycle_index = seen_states[(A, B, C)]
    non_repeating_part = coeffs[:start_of_cycle_index]
    repeating_part = coeffs[start_of_cycle_index:]
    
    return non_repeating_part, repeating_part

def main():
    """
    Main function to compute and print the continued fraction for m_{4/7}.
    """
    p, q = 4, 7
    
    # Calculate the initial number for clarity in the output
    D_val = 9 * q**2 - 4 * p**2
    num_str = f"({3*q} - sqrt({D_val})) / {2*p}"

    non_repeating, repeating = get_continued_fraction_for_quadratic_irrational(p, q)
    
    # Format the output string
    non_repeating_str = ", ".join(map(str, non_repeating[1:]))
    repeating_str = ", ".join(map(str, repeating))
    
    # The problem asks to output each number in the final equation.
    # We interpret this as clearly stating the coefficients of the continued fraction.
    print(f"The quadratic irrational number associated with m_{p}/{q} is {num_str}.")
    
    if not non_repeating: # Purely periodic
        result_str = f"[{repeating[0]}; ({', '.join(map(str, repeating[1:]))})]"
    elif len(non_repeating) == 1: # No pre-period after a0
        result_str = f"[{non_repeating[0]}; ({repeating_str})]"
    else:
        result_str = f"[{non_repeating[0]}; {non_repeating_str}, ({repeating_str})]"

    print(f"Its continued fraction is: {result_str}")
    
    # For the final answer format
    global final_answer
    final_answer = result_str


if __name__ == "__main__":
    final_answer = ""
    main()
    # The final answer format requested by the prompt
    # print(f"<<<{final_answer}>>>")