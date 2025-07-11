import math

def encode_set_to_real(S):
    """Encodes a set of natural numbers S into a real number a_S."""
    a_S = 0.0
    equation_parts = []
    for k in S:
        term = 2**(-k)
        a_S += term
        equation_parts.append(f"2**(-{k})")
    
    equation_str = " + ".join(equation_parts)
    print(f"The real number a_S for S = {S} is constructed as:\na_S = {equation_str}\na_S = {a_S}\n")
    return a_S

def check_membership_from_real(k, a_S):
    """
    Checks if k is in the set S encoded by a_S.
    This simulates the logic: k is in S iff floor(2**k * a_S) is odd.
    """
    if not isinstance(k, int) or k <= 0:
        return False, "k must be a positive integer"
        
    val = (2**k) * a_S
    floor_val = math.floor(val)
    is_odd = (floor_val % 2 == 1)
    
    # Building the 'equation' string for the output
    check_eq = f"Is floor(2**{k} * {a_S:.6f}) odd?  ->  floor({val:.6f}) = {floor_val}. Is {floor_val} odd?"
    print(check_eq, end=" ")
    return is_odd

def main():
    """
    Main function to demonstrate the encoding and decoding principle.
    We will use the set of the first 5 prime numbers as an example.
    """
    # Any subset of N can be chosen. For this demonstration, let's use a simple one.
    S = {2, 3, 5, 7, 11}
    
    # 1. Encode the set into a real number parameter a_S.
    a_S = encode_set_to_real(S)

    # 2. For any integer k, check if it's in S by 'decoding' a_S.
    print("Now, let's check for membership of various numbers k:")
    for k in range(1, 13):
        is_member = check_membership_from_real(k, a_S)
        print(f"-> {is_member}. So, k={k} is in S: {k in S}")
        assert is_member == (k in S)

if __name__ == "__main__":
    main()
