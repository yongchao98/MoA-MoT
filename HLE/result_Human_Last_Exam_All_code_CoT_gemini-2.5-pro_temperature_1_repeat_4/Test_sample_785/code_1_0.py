import numpy as np

def count_partitions(target, coeffs):
    """
    Generic function to count partitions using polynomial multiplication.
    target: The integer target sum.
    coeffs: A dictionary {dimension: num_irreps_of_that_dim}.
    """
    
    # Start with the polynomial for the constant term 1
    # which is the identity for polynomial multiplication.
    # Represents partitions of 0.
    result_poly = np.zeros(target + 1)
    result_poly[0] = 1
    
    for dim, num_irreps in coeffs.items():
        # The generating function for one irrep of dimension `dim` is 1/(1-x^dim)
        # Its polynomial representation is 1 at powers of x^dim
        irrep_poly = np.zeros(target + 1)
        for i in range(0, target + 1, dim):
            irrep_poly[i] = 1
        
        # If there are `num_irreps` of this dimension, the generating function
        # is (1/(1-x^dim))^num_irreps. We achieve this by convolving
        # the polynomial with itself `num_irreps` times.
        # Here we convolve the irrep_poly with the result_poly `num_irreps` times.
        for _ in range(num_irreps):
            result_poly = np.convolve(result_poly, irrep_poly)[:target+1]
            
    return int(round(result_poly[target]))

def main():
    # Irreducible representations of S_5
    # Dimensions: {1, 4, 5, 6}
    # Number of irreps for each dimension:
    # dim 1: 2
    # dim 4: 2
    # dim 5: 2
    # dim 6: 1
    
    target_dim = 1000
    
    # Equation:
    # n_1*1 + n_2*1 + n_3*4 + n_4*4 + n_5*5 + n_6*5 + n_7*6 = 1000
    print("The problem is equivalent to finding the number of non-negative integer solutions to the equation:")
    print("n_1*1 + n_2*1 + n_3*4 + n_4*4 + n_5*5 + n_6*5 + n_7*6 = 1000")
    
    irrep_counts = {
        1: 2,
        4: 2,
        5: 2,
        6: 1
    }

    # This is equivalent to finding the coefficient of x^1000 in the expansion of:
    # 1 / ( (1-x)^2 * (1-x^4)^2 * (1-x^5)^2 * (1-x^6)^1 )

    # A more direct way to calculate the coefficient without full polynomial multiplication
    dp = np.zeros(target_dim + 1)
    dp[0] = 1

    # For dim 1 (2 irreps)
    # This corresponds to multiplying by 1/(1-x)^2.
    # Coeffs of 1/(1-x) are all 1.
    # Coeffs of 1/(1-x)^2 are n+1.
    temp_dp = np.zeros(target_dim + 1)
    for i in range(target_dim + 1):
        for j in range(i, target_dim + 1):
            temp_dp[j] += dp[j-i]
    dp = temp_dp
    temp_dp = np.zeros(target_dim + 1)
    for i in range(target_dim + 1):
        for j in range(i, target_dim + 1):
            temp_dp[j] += dp[j-i]
    dp = temp_dp

    # For dim 4 (2 irreps)
    for _ in range(2):
        for i in range(4, target_dim + 1):
            dp[i] += dp[i - 4]
            
    # For dim 5 (2 irreps)
    for _ in range(2):
        for i in range(5, target_dim + 1):
            dp[i] += dp[i - 5]
    
    # For dim 6 (1 irrep)
    for i in range(6, target_dim + 1):
        dp[i] += dp[i - 6]

    num_orbits = int(round(dp[target_dim]))
    print(f"\nThe number of solutions, which is the number of orbits, is:")
    print(num_orbits)
    
if __name__ == "__main__":
    main()
