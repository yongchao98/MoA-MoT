def solve_ur():
    """
    Calculates and prints the minimal order u_r(n) of the Picard-Fuchs
    differential equation for the Hamiltonian V(q) = 1/2 * (q^2 - q^n)
    for n = 3, 4, ..., 12.
    """

    def get_ur(n):
        """
        Calculates u_r(n) based on the parity of n.
        The formula is n-1 for odd n and n/2 for even n.
        """
        if n % 2 == 0:
            # n is even
            return n // 2
        else:
            # n is odd
            return n - 1

    # We need to find {u_r(3), u_r(4), ..., u_r(12)}
    results = []
    for n in range(3, 13):
        ur_n = get_ur(n)
        results.append(ur_n)
        print(f"u_r({n}) = {ur_n}")
        
    # The problem asks for the set of values. We can print it as a list.
    print("\nThe set {u_r(3), u_r(4), ..., u_r(12)} is:")
    print(results)


solve_ur()