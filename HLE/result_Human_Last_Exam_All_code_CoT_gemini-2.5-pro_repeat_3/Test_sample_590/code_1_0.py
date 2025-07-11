import math

def get_degeneracy(n, k):
    """
    Calculates the degeneracy of the k-th spherical harmonic on S^(n-1).
    """
    if k == 0:
        return 1
    # We use math.comb for combinations (nCr)
    # The formula is d_k = C(n+k-1, k) - C(n+k-3, k-2)
    # C(a, b) is 0 if b < 0 or b > a.
    
    # C(n+k-1, k)
    term1 = math.comb(n + k - 1, k)
    
    # C(n+k-3, k-2)
    if k >= 2:
        term2 = math.comb(n + k - 3, k - 2)
    else:
        term2 = 0
        
    return term1 - term2

def count_positive_eigenvalues(n):
    """
    Calculates the number of positive eigenvalues for the stability operator L.
    
    Args:
        n (int): The dimension parameter from the problem description.
                 It corresponds to the sphere S^(n-1). Must be >= 2.
    """
    if not isinstance(n, int) or n < 2:
        print("Error: n must be an integer greater than or equal to 2.")
        return

    total_positive_eigenvalues = 0
    k = 0
    contributions = []
    
    # We loop through k until the condition for positive potential is no longer met.
    while True:
        # Condition for the potential V_k(rho) to be positive somewhere
        condition = k*k + (n-2)*k - n*(n-1)
        
        if condition < 0:
            degeneracy = get_degeneracy(n, k)
            contributions.append(str(degeneracy))
            total_positive_eigenvalues += degeneracy
            k += 1
        else:
            # For this k and all larger k, the potential is non-positive.
            break
            
    # As per the instruction, printing the final equation
    if not contributions:
        print("There are no positive eigenvalues found for n =", n)
        print("Final Answer: 0")
    else:
        equation_str = " + ".join(contributions)
        print(f"For n = {n}, the number of positive eigenvalues is the sum of degeneracies for k values where the potential can be positive.")
        print(f"The calculation is: {equation_str} = {total_positive_eigenvalues}")
        print(f"Final Answer: {total_positive_eigenvalues}")


# You can set the value of n here.
# As an example, we use n=3, which corresponds to a catenoid in R^4.
# The standard catenoid is in R^3, which would be n=2.
# Let's use n=3.
n = 3
count_positive_eigenvalues(n)
