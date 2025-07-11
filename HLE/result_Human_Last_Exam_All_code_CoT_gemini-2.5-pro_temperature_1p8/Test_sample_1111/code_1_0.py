import math

def demonstrate_series_behavior(p, max_terms, description):
    """
    Calculates and prints partial sums of a p-series to demonstrate its
    convergence or divergence. The equation is Sum = 1/1^p + 1/2^p + ...
    """
    print(f"--- Analyzing for {description} (p = {p}) ---")
    print(f"This corresponds to the sum: 1/1^{p} + 1/2^{p} + 1/3^{p} + ...")
    
    partial_sum = 0.0
    # Define points at which to print the running sum
    milestones = {10, 100, 1000, 10000, 100000, 1000000, 10000000}
    
    for n in range(1, max_terms + 1):
        term = 1.0 / (n**p)
        partial_sum += term
        if n in milestones:
            # Print the equation up to the current term's value, and the sum
            print(f"Sum up to n={n:<8}: {partial_sum:.6f}")
            
    print("-" * 50)
    if p > 1:
        print(f"Result: The sum converges to a finite value as n grows.")
        # For p=1.5, the sum converges to the Riemann zeta function value Zeta(1.5)
        # which is approximately 2.612375
        print(f"The infinite sum for p={p} is finite, ensuring E[T] < infinity.")
    else:
        print(f"Result: The sum appears to grow indefinitely (diverges).")
        print(f"The infinite sum for p={p} is infinite, implying E[T] = infinity.")
    print("\n")


# For k=2, the survival probability behaves like t^-1. This corresponds to a p-series with p=1.
demonstrate_series_behavior(p=1.0, max_terms=1000000, description="k=2 particles")

# For k=3, the survival probability behaves like t^(-3/2). This corresponds to a p-series with p=1.5.
demonstrate_series_behavior(p=1.5, max_terms=1000000, description="k=3 particles")