import math

def calculate_rates(T):
    """
    Calculates and prints different potential convergence rates for a given T.
    
    The problem states that in the regime T = O(e^D), the optimal rate of convergence
    for stochastic logistic regression is sought. This falls into a class of "moderately ill-posed"
    problems where the rate is faster than the general O(1/sqrt(T)) but slower than the
    O(1/T) seen in strongly convex problems. The theoretical optimal rate in this regime is
    Theta(1/T^(2/3)). This script illustrates the numerical values of these different rates.
    """
    if T <= 0:
        print("T must be a positive number.")
        return

    # Rate from option A
    rate_A = 1 / T
    
    # Rate from option B
    rate_B = 1 / (T**(2/3))
    
    # Rate from option C
    rate_C = 1 / (T**(1/2))
    
    print(f"For T = {T:.0e} samples:")
    print(f"The value for a rate of 1/T is: {rate_A:.2e}")
    print(f"The value for a rate of 1/T^(2/3) is: {rate_B:.2e}")
    print(f"The value for a rate of 1/T^(1/2) is: {rate_C:.2e}")
    print("\nAs T increases, a faster rate of convergence corresponds to a smaller value.")

# Let's demonstrate with a large number of samples.
num_samples = 1e6
calculate_rates(num_samples)
