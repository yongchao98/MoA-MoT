import math

def calculate_max_order():
    """
    Calculates the maximum order among all elements of the Smith group of S.
    """

    # Using math.comb for combinations is more efficient and direct.
    def C(n, k):
        if k < 0 or k > n:
            return 0
        return math.comb(n, k)

    # Helper function to get the prime factorization of a number.
    def get_prime_factorization(n):
        factors = {}
        d = 2
        temp_n = n
        while d * d <= temp_n:
            while (temp_n % d) == 0:
                factors[d] = factors.get(d, 0) + 1
                temp_n //= d
            d += 1
        if temp_n > 1:
            factors[temp_n] = factors.get(temp_n, 0) + 1
        return factors

    # Step 1: Define graph parameters
    n = 50  # Total elements in the set
    k = 5   # Size of the subsets (vertices)
    # Adjacency is for intersection size 3. In the Johnson scheme J(n,k),
    # the relation for intersection size k-d is denoted A_d.
    # So, 3 = k - d = 5 - d, which means d = 2.
    d = 2

    print("Step 1: Graph parameters from the Johnson scheme J(n, k)")
    print(f"n = {n}, k = {k}, adjacency intersection size = {k-d}")
    print("-" * 30)

    # Step 2: Calculate the distinct eigenvalues of the adjacency matrix A (A_2)
    print("Step 2: Distinct eigenvalues of the adjacency matrix A")
    thetas = []
    for i in range(k + 1):
        # Formula for eigenvalues of A_d in J(n,k) using Eberlein polynomials
        theta_i = sum((-1)**j * C(i, j) * C(k - i, d - j) * C(n - k - i, d - j) for j in range(d + 1))
        thetas.append(theta_i)
        print(f"theta_{i} = {theta_i}")
    print("-" * 30)

    # Step 3: Calculate the distinct eigenvalues of the Seidel matrix S = J - I - 2A
    print("Step 3: Distinct eigenvalues of the Seidel matrix S")
    v = C(n, k)
    eigenvalues = []

    # The eigenvalue corresponding to the all-ones vector
    e_0 = v - 1 - 2 * thetas[0]
    eigenvalues.append(e_0)
    print(f"e_0 = C({n}, {k}) - 1 - 2*theta_0 = {v} - 1 - 2*{thetas[0]} = {e_0}")

    # The other eigenvalues
    for i in range(1, k + 1):
        e_i = -1 - 2 * thetas[i]
        eigenvalues.append(e_i)
        print(f"e_{i} = -1 - 2*theta_{i} = -1 - 2*{thetas[i]} = {e_i}")
    print("-" * 30)

    # Step 4: The maximum order is the LCM of the absolute values of the eigenvalues.
    print("Step 4: Prime factorization of the absolute values of eigenvalues")
    abs_eigenvalues = [abs(e) for e in eigenvalues]
    
    lcm_factors = {}
    for i, num in enumerate(abs_eigenvalues):
        factors = get_prime_factorization(num)
        factor_str = ' * '.join([f'{p}^{e}' for p, e in sorted(factors.items())])
        print(f"|e_{i}| = {num} = {factor_str}")
        for p, e in factors.items():
            lcm_factors[p] = max(lcm_factors.get(p, 0), e)
    print("-" * 30)

    # Step 5: Compute the LCM from the collected prime powers
    print("Step 5: Compute the final result for the maximum order (LCM)")
    
    sorted_primes = sorted(lcm_factors.keys())
    
    terms = []
    result = 1
    for p in sorted_primes:
        power = lcm_factors[p]
        term = p**power
        terms.append(term)
        result *= term

    equation_str = " * ".join(map(str, terms))
    print("The final equation for the maximum order is:")
    print(f"{equation_str} = {result}")
    
    return result

calculate_max_order()