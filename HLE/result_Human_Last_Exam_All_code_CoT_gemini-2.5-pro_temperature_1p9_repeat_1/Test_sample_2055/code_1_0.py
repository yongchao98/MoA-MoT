import math

def prime_factorize(n):
    """
    Returns the prime factorization of n as a dictionary {prime: exponent}.
    """
    n = abs(n)
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
    Solves the problem of finding the maximum order in the Smith group of S.
    """
    # Step 1: Define parameters of the Johnson graph
    v = 50
    k = 5
    # The adjacency is defined by intersection size of 3.
    # In the language of association schemes, this is the graph A_e
    # where the intersection size is k-e, so 3 = 5-e -> e=2.
    e = 2

    # Step 2: Calculate eigenvalues of the adjacency matrix A
    # The eigenvalues mu_j of A are given by Eberlein polynomials E_e(j)
    mu = []
    print("Eigenvalues of the adjacency matrix A:")
    for j in range(k + 1):
        # mu_j = sum_{i=0 to e} (-1)^i * C(j,i) * C(k-j, e-i) * C(v-k-j, e-i)
        term1 = math.comb(j, 0) * math.comb(k - j, e - 0) * math.comb(v - k - j, e - 0) if (e-0 >= 0 and k-j >= e-0 and v-k-j >= e-0) else 0
        term2 = - (math.comb(j, 1) * math.comb(k - j, e - 1) * math.comb(v - k - j, e - 1) if (j>=1 and e-1 >= 0 and k-j >= e-1 and v-k-j >= e-1) else 0)
        term3 = math.comb(j, 2) * math.comb(k - j, e - 2) * math.comb(v - k - j, e - 2) if (j>=2 and e-2 >= 0 and k-j >= e-2 and v-k-j >= e-2) else 0
        mu_j = term1 + term2 + term3
        mu.append(mu_j)
        print(f"mu_{j} = {mu_j}")

    # Step 3: Calculate eigenvalues of the Seidel matrix S
    # S = J - I - 2A
    lambdas = []
    N = math.comb(v, k)
    
    # Eigenvalue for j=0
    lambda_0 = N - 1 - 2 * mu[0]
    lambdas.append(lambda_0)
    
    # Eigenvalues for j>0
    for j in range(1, k + 1):
        lambda_j = -1 - 2 * mu[j]
        lambdas.append(lambda_j)

    print("\nEigenvalues of the Seidel matrix S:")
    for j, l in enumerate(lambdas):
        print(f"lambda_{j} = {l}")
        
    abs_lambdas = [abs(l) for l in lambdas]

    # Step 4: Compute the LCM of the absolute values of the eigenvalues.
    print("\nTo find the maximum order, we compute the least common multiple (LCM) of the absolute values of these eigenvalues.")
    
    all_factors = {}
    for i, val in enumerate(abs_lambdas):
        factors = prime_factorize(val)
        for p, exp in factors.items():
            all_factors[p] = max(all_factors.get(p, 0), exp)

    # Compute LCM from the highest powers of primes
    lcm = 1
    sorted_primes = sorted(all_factors.keys())
    
    equation_parts = []
    for p in sorted_primes:
        exp = all_factors[p]
        term = p**exp
        lcm *= term
        equation_parts.append(f"{p}^{exp}")

    final_equation = " * ".join(equation_parts)
    print(f"\nThe maximum order is LCM(|{abs_lambdas[0]}|, |{abs_lambdas[1]}|, |{abs_lambdas[2]}|, |{abs_lambdas[3]}|, |{abs_lambdas[4]}|, |{abs_lambdas[5]}|).")
    print(f"The prime factorization of the LCM is: {final_equation}")
    print(f"The maximum order is {lcm}.")

if __name__ == "__main__":
    main()