def solve():
    """
    Calculates the number of allowed ordered pairs (a,b) with 1 <= a, b <= 1000.
    An ordered pair (a,b) is allowed if and only if gcd(a,b) = 1.
    This function counts such pairs.
    """

    # The upper limit for a and b
    limit = 1000

    # Step 1: Compute phi(k) for k from 1 to the limit using a sieve.
    # phi[k] will store Euler's totient function of k.
    phi = list(range(limit + 2))
    phi[1] = 1

    for i in range(2, limit + 1):
        # If phi[i] == i, then i is a prime number.
        if phi[i] == i:
            # Iterate through all multiples of i.
            for j in range(i, limit + 1, i):
                # Update phi[j] based on the prime factor i.
                # phi(j) = j * product_{p|j}(1 - 1/p)
                # This operation is equivalent to phi[j] = phi[j] * (i-1) / i
                phi[j] -= phi[j] // i

    # Step 2: Sum the phi values from 1 to the limit.
    # This sum gives the number of coprime pairs (a,b) with 1 <= b <= a <= limit.
    sum_phi = 0
    for i in range(1, limit + 1):
        sum_phi += phi[i]

    # Step 3: Calculate the total number of coprime pairs (a,b) where 1 <= a,b <= limit.
    # The formula is 2 * sum_phi - 1, because pairs with a!=b are counted twice
    # by symmetry, and the pair (1,1) is counted once.
    count = 2 * sum_phi - 1
    
    # As requested, output the numbers in the final equation.
    print(f"2 * {sum_phi} - 1 = {count}")

solve()