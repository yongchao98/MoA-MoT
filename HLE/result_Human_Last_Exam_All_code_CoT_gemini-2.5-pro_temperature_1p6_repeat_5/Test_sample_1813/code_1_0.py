import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def simplify_fraction(num, num_sqrt, den, D):
    """Simplifies the fractional representation by dividing by the GCD."""
    if den < 0:
        num = -num
        num_sqrt = -num_sqrt
        den = -den
    common_divisor = gcd(gcd(num, num_sqrt), den)
    return num // common_divisor, num_sqrt // common_divisor, den // common_divisor

def compute_generalized_markov_cf(p, q):
    """
    Computes the continued fraction associated with the generalized Markov number m_{p/q}.
    This corresponds to the continued fraction of alpha = (3q - 2p + sqrt(9q^2 - 4)) / (2q).
    """
    D = 9 * q**2 - 4
    
    # Initial state for alpha
    num = 3 * q - 2 * p
    num_sqrt = 1
    den = 2 * q

    # Store computed coefficients and visited states to detect cycles
    coeffs = []
    visited_states = {}
    
    # The initial state, simplified
    state = simplify_fraction(num, num_sqrt, den, D)

    # Loop until a state repeats
    while state not in visited_states:
        visited_states[state] = len(coeffs)
        
        num, num_sqrt, den = state
        
        # Calculate the integer part a_k
        # a_k = floor((num + num_sqrt * sqrt(D)) / den)
        # Using floating point for floor is safe here
        a_k = math.floor((num + num_sqrt * math.sqrt(D)) / den)
        coeffs.append(a_k)
        
        # Calculate the next state for alpha_{k+1} = 1 / (alpha_k - a_k)
        # alpha_k - a_k = (num + num_sqrt*sqrt(D))/den - a_k
        #               = (num - a_k*den + num_sqrt*sqrt(D)) / den
        
        # 1 / (alpha_k - a_k) = den / (num - a_k*den + num_sqrt*sqrt(D))
        # Rationalize the denominator:
        # Multiply by (num - a_k*den - num_sqrt*sqrt(D))
        
        next_num_part = num - a_k * den
        
        next_den = next_num_part**2 - (num_sqrt**2) * D
        next_num = den * next_num_part
        next_num_sqrt = -den * num_sqrt

        state = simplify_fraction(next_num, next_num_sqrt, next_den, D)

    # A cycle is detected. The continued fraction is periodic.
    # The period starts at the index stored in visited_states[state].
    pre_period_len = visited_states[state]
    pre_period = coeffs[:pre_period_len]
    periodic_part = coeffs[pre_period_len:]
    
    # According to theory, for a reduced quadratic irrational, the pre_period is empty.
    # Let's verify our alpha is reduced:
    # alpha = (13 + sqrt(437))/14 ~= 2.42 > 1
    # alpha' = (13 - sqrt(437))/14 ~= -0.56, so -1 < alpha' < 0
    # Since alpha is reduced, its continued fraction must be purely periodic.
    # This means pre_period should be empty.

    print("The generalized Markov number m_{4/7} is associated with the continued fraction of alpha = (13 + sqrt(437))/14.")
    print("The continued fraction is purely periodic.")
    
    final_equation = f"m_{{4/7}} -> [\\overline{{{', '.join(map(str, periodic_part))}}}]"
    print("The final equation is:")
    print(final_equation)
    
    # Per instructions, output each number in the final equation
    print("\nThe numbers in the periodic part of the continued fraction are:")
    for number in periodic_part:
        print(number)


compute_generalized_markov_cf(4, 7)