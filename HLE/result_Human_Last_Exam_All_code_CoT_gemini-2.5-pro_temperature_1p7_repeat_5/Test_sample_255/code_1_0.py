import math

def euler_totient(n):
    """
    Calculates Euler's totient function phi(n).
    For a prime power p^k, phi(p^k) = p^k - p^(k-1).
    """
    if n <= 0:
        return 0
    # The numbers we need are powers of 2.
    if (n > 0) and ((n & (n - 1)) == 0): # Check if n is a power of 2
        if n == 1:
            return 1
        return n - n // 2
    
    # General implementation for completeness, though not strictly needed here.
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

# The dimension of H^2(G,M) is given by the sum phi(d) for d in {2, 4, 8}.
d_values = [2, 4, 8]
phi_values = [euler_totient(d) for d in d_values]

# Output the explanation and the calculation
print("The dimension of the cohomology group H^2(G,M) is the sum of Euler's totient function phi(d) for d=2, 4, and 8.")
print("-" * 20)
for d, p in zip(d_values, phi_values):
    print(f"phi({d}) = {p}")

total_dimension = sum(phi_values)

# Display the final equation
equation_parts = [str(p) for p in phi_values]
print("-" * 20)
print(f"Total dimension = {' + '.join(equation_parts)} = {total_dimension}")