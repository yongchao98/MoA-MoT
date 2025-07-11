import math

def is_odd_prime(n):
    """
    Checks if a number is an odd prime.
    """
    if n <= 1 or n % 2 == 0:
        return False
    if n == 2: # This case is excluded by the previous check, but good for clarity
        return False 
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        if n % i == 0:
            return False
    return True

def find_nonabelian_filled_groups(q, m):
    """
    Finds and prints nonabelian filled groups of order 2q^m for an odd prime q.
    The classification is based on the corrected Theorem A from "A correction to
    ‘Filled groups of order pq^n’" by Glasby, Murray, and Praeger (2014).
    """
    if not isinstance(q, int) or not isinstance(m, int) or m < 1:
        print("Input error: q must be an integer and m must be a natural number (m >= 1).")
        return

    if not is_odd_prime(q):
        print(f"Input error: q must be an odd prime, but received q={q}.")
        return

    group_order = f"2*{q}^{m}"
    if m == 1:
        group_order = f"2*{q}"
    
    print(f"Searching for nonabelian filled groups of order {group_order} = {2 * (q**m)}...")

    found_groups = []
    
    # Case m=1:
    if m == 1:
        # The group is D_{2q} (Dihedral group of order 2q)
        found_groups.append(f"D_{{2*{q}}}")

    # Case m=2:
    elif m == 2:
        # One family is C_q x D_{2q}
        found_groups.append(f"C_{q} x D_{{2*{q}}}")
        
        # A special case exists for q=3
        if q == 3:
            # The group is D_18 (Dihedral group of order 18)
            found_groups.append(f"D_{{2*3^2}}")
            
    # For m > 2, the classification does not provide any such nonabelian filled groups.
    
    if found_groups:
        print("Found the following nonabelian filled group(s):")
        for group in found_groups:
            print(group)
    else:
        print("According to the known classification, no nonabelian filled groups exist for this order.")


if __name__ == '__main__':
    # --- Example Usage ---
    
    # Example 1: q=5, m=1. Order = 2 * 5^1 = 10
    print("--- Example 1: q=5, m=1 ---")
    find_nonabelian_filled_groups(q=5, m=1)
    print("\n")
    
    # Example 2: q=3, m=2. Order = 2 * 3^2 = 18
    print("--- Example 2: q=3, m=2 ---")
    find_nonabelian_filled_groups(q=3, m=2)
    print("\n")

    # Example 3: q=7, m=2. Order = 2 * 7^2 = 98
    print("--- Example 3: q=7, m=2 ---")
    find_nonabelian_filled_groups(q=7, m=2)
    print("\n")
    
    # Example 4: q=3, m=3. Order = 2 * 3^3 = 54
    print("--- Example 4: q=3, m=3 ---")
    find_nonabelian_filled_groups(q=3, m=3)
