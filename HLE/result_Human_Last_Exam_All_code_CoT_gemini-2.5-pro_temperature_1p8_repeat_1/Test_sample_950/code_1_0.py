import math

def calculate_torsion_rank(k, n):
    """
    Calculates the rank of the torsion subgroup of the integral cohomology
    of the real Grassmannian Gr_k(R^n).
    """
    
    # Let k_prime = min(k, n-k) and m_prime = max(k, n-k).
    # The rank is the number of pairs (i, j) with 1 <= j < i <= m_prime
    # such that (k_prime - i + j) is odd.
    
    k_prime = min(k, n - k)
    m_prime = max(k, n - k)
    
    print(f"For the Grassmannian Gr({k}, R^{n}):")
    print(f"We use the symmetric formula where k' = min({k}, {n-k}) = {k_prime} and m' = max({k}, {n-k}) = {m_prime}.")
    print("The rank is the number of pairs (i, j) such that 1 <= j < i <= m' and (k' - i + j) is odd.")
    print("-" * 40)
    
    count = 0
    # Iterate through all pairs (i, j) where 1 <= j < i <= m_prime
    for i in range(2, m_prime + 1):
        for j in range(1, i):
            value = k_prime - i + j
            
            print(f"Checking pair (i={i}, j={j}):")
            # This is the "equation" for each step, as requested.
            print(f"  Value = {k_prime} - {i} + {j} = {value}")
            
            if value % 2 != 0:
                print(f"  Result: {value} is odd. This pair contributes to the rank.")
                count += 1
            else:
                print(f"  Result: {value} is even. This pair does not contribute.")
            print()
            
    print("-" * 40)
    print(f"The total count of such pairs is {count}.")
    print(f"Therefore, the rank of the torsion subgroup of the integral cohomology ring is {count}.")

# Parameters for the space of 3-subspaces of R^5
k_val = 3
n_val = 5

calculate_torsion_rank(k_val, n_val)