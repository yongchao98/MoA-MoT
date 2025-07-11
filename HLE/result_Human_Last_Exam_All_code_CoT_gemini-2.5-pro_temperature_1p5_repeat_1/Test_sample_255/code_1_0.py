import math

def phi(n):
    """
    Computes Euler's totient function phi(n).
    """
    if n == 1:
        return 1
    result = n
    p = 2
    temp_n = n
    while p * p <= temp_n:
        if temp_n % p == 0:
            while temp_n % p == 0:
                temp_n //= p
            result -= result // p
        p += 1
    if temp_n > 1:
        result -= result // temp_n
    return result

# The dimension of the cohomology group H^2(G,M) is given by the sum of
# Euler's totient function for the divisors of 8 (excluding 1).
d_values = [2, 4, 8]
phi_values = [phi(d) for d in d_values]
dimension = sum(phi_values)

print("The dimension of the cohomology group H^2(G, M) is calculated by the formula:")
# Create the equation string with each number
equation_str = "dim H^2(G,M) = " + " + ".join([f"phi({d})" for d in d_values])
print(equation_str)
equation_str_values = "             = " + " + ".join([str(p) for p in phi_values])
print(equation_str_values)
print(f"             = {dimension}")
