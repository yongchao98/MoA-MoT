import math

def find_nonabelian_filled_groups_of_specific_order():
    """
    Finds which of the known nonabelian filled groups have an order of the form
    2 * q^m, where q is an odd prime and m is a natural number (m >= 1).
    """

    # Step 1: Use the known classification of nonabelian filled groups.
    # Based on the paper "A classification of finite filled groups" (Elsholtz, Lev, Shkredov, 2022),
    # there are only three such groups.
    nonabelian_filled_groups = {
        "D_10": {"order": 10},               # Dihedral group of order 10
        "C4_rtimes_C4": {"order": 16},       # A semidirect product of C4 by C4
        "(C4_x_C2)_rtimes_C2": {"order": 16} # A semidirect product of (C4 x C2) by C2
    }

    print("Identifying nonabelian filled groups of order 2*q^m (q is an odd prime, m >= 1).\n")
    print("According to current classification theorems, there are only three nonabelian filled groups.")
    print("We will check the order of each one.\n")
    
    # Step 2: Define a helper function to check the order condition.
    def check_order_form(n):
        """
        Checks if an integer n can be expressed as 2 * q^m.
        Returns a tuple (q, m) if it matches, otherwise returns None.
        """
        if n % 2 != 0 or n <= 2:
            return None
        
        k = n // 2  # k must be q^m

        if k == 1:  # m must be >= 1
            return None

        # Find all unique prime factors of k
        temp_k = k
        prime_factors = set()
        d = 2
        while d * d <= temp_k:
            if temp_k % d == 0:
                prime_factors.add(d)
                while temp_k % d == 0:
                    temp_k //= d
            d += 1
        if temp_k > 1:
            prime_factors.add(temp_k)

        # If k is a power of a prime, it will have only one prime factor.
        if len(prime_factors) != 1:
            return None

        q = prime_factors.pop()

        # The prime q must be odd.
        if q % 2 == 0:
            return None

        # Calculate m
        m = round(math.log(k, q))
        
        # Double check if k is a perfect power of q
        if q**m == k and m >= 1:
            return q, m
        else:
            return None

    # Step 3: Iterate through the groups and check them.
    found_groups = []
    for name, properties in nonabelian_filled_groups.items():
        order = properties["order"]
        print(f"--- Checking Group: {name}, Order: {order} ---")
        
        result = check_order_form(order)
        
        if result:
            q, m = result
            print(f"SUCCESS: The order {order} matches the form 2 * q^m.")
            # Final equation output as requested
            print(f"The equation is: {order} = 2 * {q}^{m}")
            found_groups.append(name)
        else:
            print(f"FAIL: The order {order} cannot be written as 2 * q^m for an odd prime q.")
        print("-" * (26 + len(name)))

    # Step 4: Print the final conclusion.
    print("\nConclusion:")
    if found_groups:
        print("The following nonabelian filled group(s) have an order of the specified form:")
        for group_name in found_groups:
            print(f"  - {group_name}")
    else:
        print("No nonabelian filled groups match the specified order form.")

if __name__ == '__main__':
    find_nonabelian_filled_groups_of_specific_order()