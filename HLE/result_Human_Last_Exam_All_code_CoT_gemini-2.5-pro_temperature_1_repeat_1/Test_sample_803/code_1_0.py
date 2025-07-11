def find_nonabelian_filled_groups(q, m):
    """
    Determines the nonabelian filled groups of order 2*q^m based on the
    known classification of such groups.

    Args:
        q (int): An odd prime number.
        m (int): A natural number (m >= 1).
    """
    # Basic input validation
    if not isinstance(q, int) or not isinstance(m, int) or m < 1:
        print("Error: q must be an integer and m must be a natural number (>= 1).")
        return
    # A simple primality test for odd numbers
    is_prime = (q > 1) and (q % 2 != 0)
    if q == 2:
      is_prime = True
    if q > 2:
        for i in range(3, int(q**0.5) + 1, 2):
            if q % i == 0:
                is_prime = False
                break
    if q != 2 and not is_prime:
        print(f"Error: The number q={q} is not an odd prime.")
        return

    order = 2 * (q**m)
    print(f"Searching for nonabelian filled groups of order 2 * {q}^{m} = {order}")

    # Based on the classification of finite filled groups, such groups
    # only exist for m=1.
    if m == 1:
        print("Conclusion: A group with these properties exists.")
        print(f"It is the dihedral group D_{2*q}, which is nonabelian and filled.")
        if q == 3:
            print("Note: For q=3, D_6 is isomorphic to the symmetric group S_3.")
    else:
        print("Conclusion: No nonabelian filled group with these properties exists for m > 1.")

if __name__ == '__main__':
    # --- Example Cases ---
    
    # Case 1: q=3, m=1. Order = 2 * 3^1 = 6
    find_nonabelian_filled_groups(3, 1)
    print("-" * 20)

    # Case 2: q=5, m=1. Order = 2 * 5^1 = 10
    find_nonabelian_filled_groups(5, 1)
    print("-" * 20)
    
    # Case 3: q=7, m=1. Order = 2 * 7^1 = 14
    find_nonabelian_filled_groups(7, 1)
    print("-" * 20)

    # Case 4: q=3, m=2. Order = 2 * 3^2 = 18
    find_nonabelian_filled_groups(3, 2)
    print("-" * 20)
