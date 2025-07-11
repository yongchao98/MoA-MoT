import math

def get_prime_factorization(n):
    """
    Returns a dictionary of prime factors of n.
    """
    factors = {}
    num = abs(n)
    d = 2
    while d * d <= num:
        while num % d == 0:
            factors[d] = factors.get(d, 0) + 1
            num //= d
        d += 1
    if num > 1:
        factors[num] = factors.get(num, 0) + 1
    return factors

def main():
    """
    Calculates the maximum order of the Smith group for the given Seidel matrix.
    """
    v = 50
    k = 5

    # Use math.comb for combinations
    comb = math.comb

    # Step 1: Calculate the eigenvalues of the adjacency matrix A
    print(f"The graph is defined on the {k}-subsets of a {v}-set.")
    print("Adjacency is defined by intersection size 3.")
    print("\nStep 1: Calculating the eigenvalues of the adjacency matrix A...")
    
    lambdas = []
    for j in range(k + 1):
        # Formula for eigenvalues of the graph A_i in Johnson scheme J(v,k)
        # Here i=2, since intersection size is k-2 = 3.
        val = sum(
            ((-1)**t) * comb(j, t) * comb(k - j, 2 - t) * comb(v - k - j, 2 - t)
            for t in range(3)
        )
        lambdas.append(val)
        print(f"  lambda_{j} = {val}")

    # Step 2: Calculate the eigenvalues of the Seidel matrix S
    print("\nStep 2: Calculating the eigenvalues of the Seidel matrix S...")
    n = comb(v, k)
    
    thetas = []
    # The eigenvalue corresponding to the all-ones vector
    theta_0 = n - 1 - 2 * lambdas[0]
    thetas.append(theta_0)
    print(f"  theta_0 = n - 1 - 2*lambda_0 = {n} - 1 - 2*{lambdas[0]} = {theta_0}")

    # The other eigenvalues
    for j in range(1, k + 1):
        theta_j = -(1 + 2 * lambdas[j])
        thetas.append(theta_j)
        print(f"  theta_{j} = -(1 + 2*lambda_{j}) = -(1 + 2*{lambdas[j]}) = {theta_j}")
        
    distinct_abs_thetas = sorted(list(set(abs(t) for t in thetas)))

    # Step 3: Find the prime factorizations of the absolute values of distinct eigenvalues
    print("\nStep 3: Finding the prime factorization of each distinct absolute eigenvalue...")
    
    all_prime_factors = {}
    for t_abs in distinct_abs_thetas:
        factors = get_prime_factorization(t_abs)
        factor_str = " * ".join([f"{p}^{e}" for p, e in factors.items()])
        print(f"  |{t_abs}| = {factor_str}")
        for p, e in factors.items():
            if p not in all_prime_factors or e > all_prime_factors[p]:
                all_prime_factors[p] = e

    # Step 4: Calculate the LCM
    print("\nStep 4: Calculating the maximum order, which is the LCM of the distinct absolute eigenvalues.")
    
    lcm = 1
    lcm_str_parts = []
    for p, e in sorted(all_prime_factors.items()):
        lcm *= (p ** e)
        if e > 1:
            lcm_str_parts.append(f"{p}^{e}")
        else:
            lcm_str_parts.append(f"{p}")
            
    lcm_eq_str = " * ".join(lcm_str_parts)
    print(f"The final calculation is:\n{lcm_eq_str}")
    
    print(f"\nThe maximum order among all elements of the Smith group of S is {lcm}.")
    
    # Final answer in the required format
    print(f"\n<<<{lcm}>>>")

if __name__ == "__main__":
    main()