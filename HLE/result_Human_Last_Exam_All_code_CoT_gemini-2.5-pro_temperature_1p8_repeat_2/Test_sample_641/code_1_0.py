import math

def solve():
    """
    Calculates the number of involutions in PSU(4, 997).
    """
    q = 997

    # Part 1: Number of elements g such that g^2=I and g is not in the center.
    # This corresponds to the conjugacy class of involutions with a 2-dimensional (-1)-eigenspace.
    # Formula: q^4 * (q^2 - q + 1) * (q^2 + 1)
    q2 = q * q
    q4 = q2 * q2
    
    n1 = q4 * (q2 - q + 1) * (q2 + 1)
    
    # Part 2: Number of elements g such that g^2 = -I.
    # For q = 997 (where q = 1 mod 4), these form one conjugacy class.
    # Formula: q^2 * (q^3 + 1)
    q3 = q2 * q
    
    n2 = q2 * (q3 + 1)
    
    # The total number of preimages in SU(4,q) is the sum of these two counts.
    total_preimages = n1 + n2
    
    # The size of the center of SU(4, q) is d = gcd(4, q+1) = gcd(4, 998) = 2.
    d = 2
    
    # The number of involutions in PSU(4,q) is the number of preimages divided by d.
    num_involutions = total_preimages // d
    
    print(f"The number of involutions is calculated using the formula: (N1 + N2) / d")
    print(f"Here, q = {q}")
    print(f"d = gcd(4, q+1) = {d}")
    print(f"N1 (elements g with g^2=I, g not in Center) = q^4(q^2-q+1)(q^2+1) = {n1}")
    print(f"N2 (elements g with g^2=-I) = q^2(q^3+1) = {n2}")
    print(f"Total number of involutions in PSU(4, {q}) = ({n1} + {n2}) / {d}")
    print(f"Final answer: {num_involutions}")

solve()