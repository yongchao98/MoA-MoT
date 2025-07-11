import math

def get_multiplicative_order(n, modulus):
    """
    Calculates the multiplicative order of n modulo modulus.
    This is the smallest k > 0 such that n^k = 1 (mod modulus).
    """
    if math.gcd(n, modulus) != 1:
        return 0
    k = 1
    power = n % modulus
    while power != 1:
        power = (power * n) % modulus
        k += 1
    return k

def calculate_scl_in_G(p, q):
    """
    Calculates the stable commutator length for an element in G
    with rational rotation number p/q.
    """
    # The scl is 0 for dyadic rationals. q must be odd.
    if q % 2 == 0:
        return 0, 1

    # Step 1: Find the length of the repeating binary block.
    # This is the multiplicative order of 2 modulo q.
    k = get_multiplicative_order(2, q)

    # Step 2: Find the repeating binary block.
    # The fraction p/q is equal to A / (2^k - 1) for some integer A.
    # A is the integer whose k-bit binary representation is the repeating block.
    numerator_A = p * (2**k - 1) // q
    binary_block = bin(numerator_A)[2:].zfill(k)

    # Step 3: Calculate the maximal fluctuation of the partial sums.
    # We treat '0' as -1 and '1' as +1.
    s = 0
    min_s = 0
    max_s = 0
    for bit in binary_block:
        if bit == '1':
            s += 1
        else:
            s -= 1
        
        if s < min_s:
            min_s = s
        if s > max_s:
            max_s = s
            
    max_fluctuation = max_s - min_s

    # Step 4: Apply the scl formula for G.
    # scl = max_fluctuation / (2 * (2^k - 1))
    scl_numerator = max_fluctuation
    scl_denominator = 2 * (2**k - 1)
    
    return scl_numerator, scl_denominator

def main():
    """
    Main function to solve the problem.
    """
    # Rotation numbers for g and h
    # tau_g = 2/27, tau_h = 16/27

    # We need scl(g^2) and scl(h^2).
    # tau(g^2) = 2 * (2/27) = 4/27
    p1, q1 = 4, 27
    
    # tau(h^2) = 2 * (16/27) = 32/27. We can use its value mod 1, which is 5/27.
    p2, q2 = 5, 27

    # Calculate scl_G(g^2)
    scl_g2_num, scl_g2_den = calculate_scl_in_G(p1, q1)

    # Calculate scl_G(h^2)
    scl_h2_num, scl_h2_den = calculate_scl_in_G(p2, q2)

    # The final scl is sqrt(scl_G(g^2) * scl_G(h^2)).
    # In this case, the two values are identical, so the result is just one of them.
    final_scl_num = scl_g2_num
    final_scl_den = scl_g2_den

    print("The stable commutator length of [g_1, h_2] is given by the equation:")
    print(f"scl = sqrt(scl_G(g^2) * scl_G(h^2))")
    print(f"scl_G(g^2) for rotation number 4/27 is: {scl_g2_num} / {scl_g2_den}")
    print(f"scl_G(h^2) for rotation number 5/27 is: {scl_h2_num} / {scl_h2_den}")
    print("Since these values are the same, the final result is:")
    print(f"{final_scl_num} / {final_scl_den}")

if __name__ == "__main__":
    main()