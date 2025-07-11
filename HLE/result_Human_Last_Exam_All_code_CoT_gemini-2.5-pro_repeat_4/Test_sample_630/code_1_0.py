import math

def compare_rates(T):
    """
    Calculates the values of different convergence rate functions for a given T.
    """
    print(f"For T = {T:.1e}:")
    
    # Rate from option A
    rate_A = 1 / T
    print(f"Rate Theta(1/T): {rate_A:.2e}")
    
    # Rate from option B
    rate_B = 1 / (T**(2/3))
    print(f"Rate Theta(1/T^(2/3)): {rate_B:.2e}")
    
    # Rate from option C
    rate_C = 1 / math.sqrt(T)
    print(f"Rate Theta(1/T^(1/2)): {rate_C:.2e}")
    
    # Rate derived from analysis
    derived_rate = math.log(T) / math.sqrt(T)
    print(f"Derived Rate Omega(log(T)/T^(1/2)): {derived_rate:.2e}")
    
    # We can see the derived rate is larger (slower) than Theta(1/sqrt(T)).
    # For example, let's check the ratio.
    ratio = derived_rate / rate_C
    print(f"\nRatio of (Derived Rate) / (Rate C) = log(T) = {ratio:.2f}")

# Execute for a large value of T
compare_rates(1_000_000_000)
