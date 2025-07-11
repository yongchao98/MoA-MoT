import math

def is_prime(n):
    """Helper function to check if a number is prime."""
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def identify_nonabelian_filled_group(q, m):
    """
    Identifies the nonabelian filled group of order 2*q^m for an odd
    prime q and natural number m, based on the known classification theorem.
    """
    # 1. Validate the inputs
    if not isinstance(q, int) or not isinstance(m, int):
        print("Error: q and m must be integers.")
        return
    if m < 1:
        print(f"Error: m must be a natural number (m >= 1). Provided m={m}.")
        return
    if q % 2 == 0 or not is_prime(q):
        print(f"Error: q must be an odd prime. Provided q={q} is not.")
        return

    # 2. Apply the classification theorem
    # The nonabelian filled groups of order 2q^m are the dihedral groups D_{2n}
    # where n = q^m.
    n = q**m
    order = 2 * n

    # 3. Print the result
    print(f"For q = {q} and m = {m}:")
    print(f"The group order is 2 * {q}^{m} = {order}.")
    print("Based on the classification of filled groups, the unique nonabelian filled group of this order is the Dihedral Group.")
    print(f"Name: D_{order} (also written as D_{n} where n is the number of rotational symmetries, n={n})")
    print("The standard presentation of this group is:")
    print(f"<r, s | r^{n} = s^2 = 1, srs = r^-1>")

# --- Example ---
# Let's identify the group for an example case: q=5 (an odd prime) and m=1 (a natural number).
q_example = 5
m_example = 1
identify_nonabelian_filled_group(q_example, m_example)
