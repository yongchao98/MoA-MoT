import math

def euler_phi(n):
    """Computes Euler's totient function phi(n)."""
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

def get_prime_factorization(n):
    """Returns a dictionary of prime factors of n."""
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

# Conductor and order
d = 53599
target_order = 6

# Step 1: Analyze the conductor d.
# d = 53599 = 7 * 13 * 19 * 31. It is a square-free product of 4 distinct primes.
# A character is primitive if its components are non-principal (order > 1).
factors = get_prime_factorization(d)
num_primes = len(factors)

# Step 2: Apply the condition on the order (lcm of component orders = 6)
# using the Principle of Inclusion-Exclusion.
# The number of characters of order k is phi(k), provided k divides phi(p) = p-1.
# All prime factors p in {7, 13, 19, 31} have p-1 divisible by 6.

# Part 1: Count tuples where each component's order divides 6 (primitive orders: 2, 3, 6).
phi_2 = euler_phi(2)
phi_3 = euler_phi(3)
phi_6 = euler_phi(6)
choices_div_6 = phi_2 + phi_3 + phi_6
total_div_6 = choices_div_6 ** num_primes

# Part 2: Count tuples where each component's order divides 3 (primitive order: 3).
choices_div_3 = phi_3
total_div_3 = choices_div_3 ** num_primes

# Part 3: Count tuples where each component's order divides 2 (primitive order: 2).
choices_div_2 = phi_2
total_div_2 = choices_div_2 ** num_primes

# Final Calculation
final_count = total_div_6 - total_div_3 - total_div_2

# Print the final equation as requested
print("The number of primitive Dirichlet characters of conductor 53599 and order 6 is calculated as follows:")
print("Let N be the number. The conductor is d = 53599 = 7 * 13 * 19 * 31, a product of 4 primes.")
print("The number of primitive components with order dividing 6 is phi(2) + phi(3) + phi(6) = {} + {} + {} = {}.".format(phi_2, phi_3, phi_6, choices_div_6))
print("The number of primitive components with order dividing 3 is phi(3) = {}.".format(phi_3))
print("The number of primitive components with order dividing 2 is phi(2) = {}.".format(phi_2))
print("By inclusion-exclusion, the total number N is:")
print("N = (choices for order dividing 6)^4 - (choices for order dividing 3)^4 - (choices for order dividing 2)^4")
print("N = ({})^{} - ({})^{} - ({})^{}".format(choices_div_6, num_primes, choices_div_3, num_primes, choices_div_2, num_primes))
print("N = {} - {} - {}".format(total_div_6, total_div_3, total_div_2))
print("N = {}".format(final_count))