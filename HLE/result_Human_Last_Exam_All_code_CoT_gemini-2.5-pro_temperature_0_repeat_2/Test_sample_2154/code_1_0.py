def solve():
    """
    This function calculates the minimal order of the Picard-Fuchs differential equation u_r(n)
    for the Hamiltonian V(q) = 1/2 * (q^2 - q^n) for n from 3 to 12.
    """

    def u_r(n):
        """
        Calculates u_r(n) based on the parity of n.
        - If n is odd, the order is n - 1.
        - If n is even, the order is n / 2 due to symmetry.
        """
        if n % 2 == 1:  # n is odd
            return n - 1
        else:  # n is even
            return n // 2

    results = []
    print("The values for the sequence {u_r(3), u_r(4), ..., u_r(12)} are:")
    for n in range(3, 13):
        order = u_r(n)
        results.append(order)
        # The problem asks to output each number in the final equation.
        # We interpret this as printing each computed value clearly.
        print(f"u_r({n}) = {order}")

    # For the final answer format, we also print the list of values.
    # print(f"\nThe complete set is: {results}")

solve()