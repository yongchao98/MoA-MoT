import math

def calculate_scl():
    """
    This function calculates the stable commutator length as requested.
    """

    # Step 1: Define the parameters for the rational rotations g and h.
    # g corresponds to rotation by p1/q1
    p1, q1 = 2, 27
    # h corresponds to rotation by p2/q2
    p2, q2 = 16, 27

    print("Step 1: Calculate the scl for g, corresponding to translation by {}/{}.".format(p1, q1))
    # Check if p1 and q1 are coprime
    if math.gcd(p1, q1) != 1:
        print("Error: p1 and q1 are not coprime.")
        return

    # Calculate the signature of the (p1, q1) torus knot
    sig1 = -(p1 - 1) * (q1 - 1)
    # Calculate the scl for g
    scl1_num = abs(sig1)
    scl1_den = 8
    scl1 = scl1_num / scl1_den
    
    print("The signature of the T({},{}) torus knot is:".format(p1, q1))
    print("  \u03C3(T{},{}) = -({} - 1) * ({} - 1) = -({}) * ({}) = {}".format(p1, q1, p1, q1, p1-1, q1-1, sig1))
    print("The stable commutator length of g is:")
    print("  scl(g) = |{}| / 8 = {}/{} = {}".format(sig1, scl1_num, scl1_den, scl1))
    print("-" * 20)

    print("Step 2: Calculate the scl for h, corresponding to translation by {}/{}.".format(p2, q2))
    # Check if p2 and q2 are coprime
    if math.gcd(p2, q2) != 1:
        print("Error: p2 and q2 are not coprime.")
        return
        
    # Calculate the signature of the (p2, q2) torus knot
    sig2 = -(p2 - 1) * (q2 - 1)
    # Calculate the scl for h
    scl2_num = abs(sig2)
    scl2_den = 8
    scl2 = scl2_num / scl2_den

    print("The signature of the T({},{}) torus knot is:".format(p2, q2))
    print("  \u03C3(T{},{}) = -({} - 1) * ({} - 1) = -({}) * ({}) = {}".format(p2, q2, p2, q2, p2-1, q2-1, sig2))
    print("The stable commutator length of h is:")
    print("  scl(h) = |{}| / 8 = {}/{} = {}".format(sig2, scl2_num, scl2_den, scl2))
    print("-" * 20)
    
    # Step 3: Compute the total stable commutator length.
    total_scl = scl1 + scl2
    total_scl_num = scl1_num * scl2_den + scl2_num * scl1_den
    total_scl_den = scl1_den * scl2_den
    # simplify the fraction
    common_divisor = math.gcd(total_scl_num, total_scl_den)
    total_scl_num_simple = total_scl_num // common_divisor
    total_scl_den_simple = total_scl_den // common_divisor

    print("Step 3: The final scl is the sum of the individual scl values.")
    print("scl(g\u2081h\u2082) = scl(g) + scl(h)")
    print("Final calculation:")
    print("  {} + {} = {}/{} + {}/{} = {}/{} = {}".format(
        scl1, scl2, scl1_num, scl1_den, scl2_num, scl2_den,
        scl1_num+scl2_num, scl1_den,
        total_scl
    ))

calculate_scl()