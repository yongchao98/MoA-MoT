import math

memo = {
    (0, 1): 1,
    (1, 1): 1,
    (1, 2): 2
}

def extended_gcd(a, b):
    """Returns (gcd, x, y) such that a*x + b*y = gcd."""
    if a == 0:
        return (b, 0, 1)
    d, x1, y1 = extended_gcd(b % a, a)
    x = y1 - (b // a) * x1
    y = x1
    return d, x, y

def find_parents(p, q):
    """
    Finds parent fractions (p1/q1, p2/q2) for p/q in the Farey tree.
    Solves the Diophantine equation q*x - p*y = 1.
    """
    if p == 1 and q == 2:
        return (0, 1), (1, 1)

    # We need to solve q*x - p*y = 1 for p1=x, q1=y
    # This is equivalent to p*y - q*x = -1
    gcd, y, x_neg, = extended_gcd(p, q) # ax+by=gcd -> py-qx=-1 requires different variables
    
    # We want to solve p*s - q*r = 1 (or -1)
    # Using EEA on (q,p): qx - py = 1
    # Note: p, q must be coprime, so gcd is 1.
    g, x, y_neg = extended_gcd(q, -p)
    
    # We have q*x + (-p)*y_neg = 1 => q*x - p*y_neg = 1
    # One solution is (r, s) = (x, y_neg)
    # We need 0 < r < p and 0 < s < q
    # General solution: r_k = x - k*p, s_k = y_neg + k*q
    # We need to find k. x - k*p > 0 => k < x/p.
    k = math.floor((x-1)/p)
    
    r1 = x - k * p
    s1 = y_neg + k * q
    
    # Check for validity
    if not (0 < r1 < p and 0 < s1 < q):
        # if qr-ps=-1, parents are (r,s) and (p-r, q-s)
        # if ps-qr=-1, parents are (p-r,q-s) and (r,s)
        # let's find a solution to ps-qr=1 -> p*s - q*r = 1
        g, s, r_neg = extended_gcd(p, -q)
        r = -r_neg
        k = math.floor((s-1)/q)
        s1 = s - k*q
        r1 = r + k*p
        
    p1, q1 = r1, s1
    p2, q2 = p - p1, q - q1

    return (p1, q1), (p2, q2)


def get_m(p, q):
    """
    Recursively computes the generalized Markov number m_{p/q} using a Farey tree recurrence.
    """
    if (p, q) in memo:
        return memo[(p, q)]

    # Ensure fraction is in simplest form and p < q
    common_divisor = math.gcd(p, q)
    p, q = p // common_divisor, q // common_divisor
    if p > q:
        p, q = q, p
    if p/q > 0.5:
        p,q = q-p, q # Use symmetry m_{p/q} = m_{1-p/q}

    if (p, q) in memo:
        return memo[(p, q)]
    
    # Find parent fractions
    (p1, q1), (p2, q2) = find_parents(p, q)

    # Find grandparent fraction
    gp, gq = abs(p1 - p2), abs(q1 - q2)

    # Recursive calls
    m1 = get_m(p1, q1)
    m2 = get_m(p2, q2)
    m0 = get_m(gp, gq)

    # The recurrence relation
    result = 3 * m1 * m2 - m0
    memo[(p, q)] = result
    return result

if __name__ == '__main__':
    p_target, q_target = 4, 7
    
    # Get values needed for the final equation
    parents = find_parents(p_target, q_target)
    (p1, q1), (p2, q2) = parents
    
    grand_p, grand_q = abs(p1 - p2), abs(q1 - q2)

    m_target = get_m(p_target, q_target)
    m1 = get_m(p1, q1)
    m2 = get_m(p2, q2)
    m0 = get_m(grand_p, grand_q)
    
    print(f"The generalized Markov number for the rational {p_target}/{q_target} is m_{{{p_target}/{q_target}}}.")
    print("This number is found using a recursive formula on the Farey tree.")
    print(f"The parents of {p_target}/{q_target} are {p1}/{q1} and {p2}/{q2}.")
    print(f"The 'grandparent' is {grand_p}/{grand_q}.")
    print("\nThe corresponding Markov numbers are:")
    print(f"m_{{{p1}/{q1}}} = {m1}")
    print(f"m_{{{p2}/{q2}}} = {m2}")
    print(f"m_{{{grand_p}/{grand_q}}} = {m0}")
    print("\nThe final computation is:")
    print(f"m_{{{p_target}/{q_target}}} = 3 * m_{{{p1}/{q1}}} * m_{{{p2}/{q2}}} - m_{{{grand_p}/{grand_q}}}")
    print(f"m_{{{p_target}/{q_target}}} = 3 * {m1} * {m2} - {m0} = {m_target}")
    
    final_markov_number = m_target
    # This checks if the result is a valid Markov number triple with its parents.
    is_valid_triple = m1**2 + m2**2 + m_target**2 == 3 * m1 * m2 * m_target
    #print(f"\nVerification: ({m1}, {m2}, {m_target}) is a Markov triple: {is_valid_triple}")
    
