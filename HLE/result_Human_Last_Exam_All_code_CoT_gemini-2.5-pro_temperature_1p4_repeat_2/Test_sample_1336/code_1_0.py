import math

def calculate_smooth_coverings(p):
    """
    Calculates the total number of smooth coverings for a given prime p > 5.

    Args:
        p: A prime number greater than 5.
    """
    if p <= 5 or not all(p % i for i in range(2, int(math.sqrt(p)) + 1)):
        print(f"Error: p must be a prime number greater than 5. Received p={p}.")
        return

    # The number of irreducible characters (or blocks in char 0) of SL(2, p) is p + 4.
    num_chars_sl2p = p + 4
    
    # The number of irreducible characters (or blocks in char 0) of PSL(2, p) is (p + 5) / 2.
    # These correspond to the non-faithful characters of SL(2, p).
    # We use integer division since p is an odd prime, so p+5 is even.
    num_chars_psl2p = (p + 5) // 2

    # The number of smooth coverings is interpreted as the number of "faithful" characters of SL(2, p),
    # which is the total number of characters of SL(2, p) minus the number of characters of PSL(2, p).
    total_coverings = num_chars_sl2p - num_chars_psl2p

    print(f"For p = {p}:")
    print(f"Step 1: The total number of irreducible characters of SL(2, p) is {p} + 4 = {num_chars_sl2p}.")
    print(f"Step 2: The number of non-faithful characters (characters of PSL(2, p)) is ({p} + 5) / 2 = {num_chars_psl2p}.")
    print(f"Step 3: The total number of smooth coverings is the difference: {num_chars_sl2p} - {num_chars_psl2p} = {total_coverings}.")
    print(f"The final formula is (p+3)/2, which for p={p} is ({p}+3)/2 = {total_coverings}.")

# Example calculation for a given prime p > 5, for example p = 11.
p = 11
calculate_smooth_coverings(p)