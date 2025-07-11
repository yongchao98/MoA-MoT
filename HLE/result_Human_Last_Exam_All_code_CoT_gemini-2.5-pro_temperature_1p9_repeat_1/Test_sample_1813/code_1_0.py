import math

def compute_markov_cf():
    """
    Computes and prints the details of the generalized Markov number m_4/7
    and its associated continued fraction.
    """
    # The generalized Markov number m_4/7 is given by the expression:
    # m = (A + sqrt(B)) / C
    A = 4
    B = 13
    C = 2

    print(f"The generalized Markov number m_4/7 is given by the equation m = (A + sqrt(B)) / C.")
    print(f"For this problem, the values are:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")
    print(f"So, m_4/7 = ({A} + sqrt({B})) / {C}")
    
    # --- Compute the continued fraction of m_4/7 ---
    # We will represent our number as (num_a + num_b*sqrt(d))/den
    d = B
    num_a = A
    num_b = 1
    den = C

    coeffs = []
    
    # Compute enough terms to see the repeating pattern. 15 is sufficient.
    for _ in range(15):
        # Calculate the integer part of the current number
        a_i = math.floor((num_a + num_b * math.sqrt(d)) / den)
        coeffs.append(a_i)
        
        # Subtract the integer part
        # xi_new = 1 / (xi - a_i)
        # 1. temp = xi - a_i = ((num_a - a_i*den) + num_b*sqrt(d))/den
        temp_num_a = num_a - a_i * den
        
        # 2. Invert and rationalize
        # 1/temp = den / (temp_num_a + num_b*sqrt(d))
        #        = den * (temp_num_a - num_b*sqrt(d)) / (temp_num_a^2 - num_b^2*d)
        new_den = temp_num_a**2 - num_b**2 * d
        new_num_a = den * temp_num_a
        new_num_b = -den * num_b
        
        # 3. Simplify the new representation by dividing by GCD
        common_divisor = math.gcd(math.gcd(new_num_a, new_num_b), new_den)
        num_a = new_num_a // common_divisor
        num_b = new_num_b // common_divisor
        den = new_den // common_divisor
        
        # 4. Conventionally, the denominator should be positive
        if den < 0:
            den = -den
            num_a = -num_a
            num_b = -num_b

    # Print the continued fraction result
    # The computed fraction is [3; 1, 4, 14, 4, 4, 14, 4, ...], which is [3; 1, overline(4, 14, 4)]
    print("\nThe continued fraction for m_4/7 is [c0; c1, c2, c3, ...]")
    print("The coefficients are:")
    
    # Create the output string like [3; 1, 4, 14, ...]
    cf_string = f"[{coeffs[0]}; {', '.join(map(str, coeffs[1:]))}, ...]"
    print(cf_string)


compute_markov_cf()