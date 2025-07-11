import math

def get_prime_factorization(n):
    """
    Computes the prime factorization of a positive integer n.
    Returns a dictionary mapping each prime factor to its exponent.
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def main():
    """
    Main function to solve the problem.
    """
    # Parameters from the problem description
    v = 50
    k = 5

    # Step 1: Compute eigenvalues of the adjacency matrix A for the graph where
    # adjacency is defined by intersection size 3. In the Johnson scheme J(v,k),
    # this corresponds to the matrix A_2.
    # The eigenvalues P_j2 are given by the formula:
    # P_j2 = C(k-j, 2)C(v-k-j, 2) - j(k-j)(v-k-j) + C(j, 2)
    P_j2 = []
    print("--- Step 1: Calculating eigenvalues of the graph's adjacency matrix A ---")
    for j in range(k + 1):
        term1 = math.comb(k - j, 2) * math.comb(v - k - j, 2) if k - j >= 2 and v - k - j >= 2 else 0
        term2 = j * (k - j) * (v - k - j) if k - j >= 0 and v - k - j >= 0 else 0
        term3 = math.comb(j, 2) if j >= 2 else 0
        
        val = term1 - term2 + term3
        P_j2.append(val)
        print(f"Eigenvalue P_{j}2 = {val}")

    # Step 2: Compute eigenvalues of the Seidel matrix S = J - I - 2A
    print("\n--- Step 2: Calculating eigenvalues of the Seidel matrix S ---")
    n = math.comb(v, k)
    theta = []
    # The eigenvalue for the all-ones eigenvector (j=0)
    theta_0 = n - 1 - 2 * P_j2[0]
    theta.append(theta_0)
    print(f"Eigenvalue theta_0 = C(50, 5) - 1 - 2 * P_02 = {n} - 1 - 2 * {P_j2[0]} = {theta_0}")

    # Eigenvalues for eigenvectors orthogonal to the all-ones vector (j>0)
    for j in range(1, k + 1):
        val = -(1 + 2 * P_j2[j])
        theta.append(val)
        print(f"Eigenvalue theta_{j} = -(1 + 2 * P_{j}2) = -(1 + 2 * {P_j2[j]}) = {val}")
    
    # Step 3: Get absolute values of eigenvalues and find their prime factorizations
    print("\n--- Step 3: Prime factorization of the absolute values of eigenvalues ---")
    abs_theta = [abs(t) for t in theta]
    factorizations = []
    for i, t in enumerate(abs_theta):
        factors = get_prime_factorization(t)
        factorizations.append(factors)
        factor_str = " * ".join([f"{p}^{e}" for p, e in sorted(factors.items())])
        print(f"|theta_{i}| = {t} = {factor_str}")

    # Step 4: Find the maximum exponent for each prime factor
    print("\n--- Step 4: Determining the maximum exponent for each prime base ---")
    max_exponents = {}
    for factors in factorizations:
        for p, e in factors.items():
            if p not in max_exponents or e > max_exponents[p]:
                max_exponents[p] = e

    sorted_primes = sorted(max_exponents.keys())
    max_exp_str_list = []
    for p in sorted_primes:
        max_exp_str_list.append(f"{p}^{max_exponents[p]}")
    print(f"The structure of the largest invariant factor d_n is the product of these prime powers:")
    print(" * ".join(max_exp_str_list))


    # Step 5: Calculate the final result, which is the largest invariant factor d_n
    print("\n--- Step 5: Calculating the final result ---")
    d_n = 1
    terms = []
    for p in sorted_primes:
        term_val = p**max_exponents[p]
        terms.append(str(term_val))
        d_n *= term_val

    print(f"The final product is d_n = {' * '.join(terms)}")
    print(f"\nThe maximum order among all elements of the Smith group of S is:")
    print(d_n)

if __name__ == "__main__":
    main()