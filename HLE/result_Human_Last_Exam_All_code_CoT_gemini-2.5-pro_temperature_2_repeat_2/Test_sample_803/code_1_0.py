def solve_for_filled_groups():
    """
    Identifies nonabelian filled groups with order of the form 2*q^m
    based on the known classification of finite filled groups.
    """

    # According to the literature, the nonabelian finite filled groups are D_10 and SmallGroup(20,3).
    # We represent them by their names and orders.
    known_nonabelian_filled_groups = {
        "The dihedral group of order 10 (D_10)": 10,
        "The semidirect product C_5 x| C_4 (SmallGroup(20,3))": 20
    }

    # This function checks if a number n can be expressed as q^m for an odd prime q.
    def get_odd_prime_power_form(n):
        """
        If n = q^m for an odd prime q and integer m >= 1, returns (q, m).
        Otherwise, returns None.
        """
        if n <= 1:
            return None

        # Special case: check if n is a power of an even prime (2)
        if n > 0 and (n & (n - 1) == 0): # Check if n is a power of 2
             return None

        # Find a prime factor
        d = 3
        temp_n = n
        while d * d <= temp_n:
            if temp_n % d == 0:
                # d is a prime factor. Check if n is a power of d.
                while temp_n % d == 0:
                    temp_n //= d
                if temp_n == 1:
                    m = 0
                    power_val = 1
                    while power_val < n:
                        power_val *= d
                        m += 1
                    return (d, m)
                else:
                    # n has multiple prime factors
                    return None
            d += 2

        # If we reach here, n must be a prime itself (since we checked for factors up to sqrt(n))
        if n > 2: # It's an odd prime
            return (n, 1)

        return None


    found_groups = []
    print("Searching for nonabelian filled groups of order 2*q^m (q is odd prime, m is natural number)...")
    print("-" * 70)

    for name, order in known_nonabelian_filled_groups.items():
        # The order must be of the form 2 * q^m. So, order/2 must be q^m.
        if order % 2 == 0:
            target = order // 2
            result = get_odd_prime_power_form(target)

            if result:
                q, m = result
                group_info = {
                    "name": name,
                    "order": order,
                    "q": q,
                    "m": m
                }
                found_groups.append(group_info)
                print(f"Found a potential match: {name}")
                print(f"  Order = {order}. Checking if it fits the form 2 * q^m.")
                print(f"  Order / 2 = {target}. Checking if {target} is a power of an odd prime.")
                print(f"  Yes, {target} = {q}^{m}. This fits the criteria.\n")
            else:
                print(f"Checking group: {name}")
                print(f"  Order = {order}. Order / 2 = {target}.")
                print(f"  {target} is not a power of a single odd prime. This group is not a solution.\n")

    print("-" * 70)
    if found_groups:
        print("Conclusion: The only nonabelian filled group of order 2*q^m is:")
        for group in found_groups:
            order = group['order']
            q = group['q']
            m = group['m']
            print(f"Group: {group['name']}")
            print(f"This is because its order, {order}, can be expressed in the required form.")
            # Final required output: each number in the final equation
            print(f"The equation for the order is: {order} = 2 * {q}**{m}")
    else:
        print("Conclusion: No nonabelian filled groups match the specified order criteria.")

if __name__ == "__main__":
    solve_for_filled_groups()