import math

def get_prime_factorization(n):
    """
    Returns a dictionary representing the prime factorization of n.
    Example: get_prime_factorization(12) returns {2: 2, 3: 1} for 12 = 2^2 * 3^1
    """
    factors = {}
    # Check for factor 2
    while n % 2 == 0:
        factors[2] = factors.get(2, 0) + 1
        n = n // 2
    # Check for odd factors
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        while n % i == 0:
            factors[i] = factors.get(i, 0) + 1
            n = n // i
    # If n is a prime number greater than 2
    if n > 2:
        factors[n] = factors.get(n, 0) + 1
    return factors

# Number of dice and sides per die
num_dice = 100
num_sides = 6

# Step 1: The size of the sample space is num_sides ^ num_dice.
# We need its prime factorization: (p1^a1 * p2^a2 * ...)^num_dice
# This is equal to p1^(a1*num_dice) * p2^(a2*num_dice) * ...

# Step 2: Get the prime factorization of the number of sides.
side_factors = get_prime_factorization(num_sides)

# Step 3: Calculate the exponents for the sample space size factorization.
space_exponents = {base: exponent * num_dice for base, exponent in side_factors.items()}

# Step 4: The maximum number of independent events is the sum of these exponents.
max_m = sum(space_exponents.values())

# Extract the numbers for the final equation output
exponents = list(space_exponents.values())
exponent_p1 = exponents[0]
exponent_p2 = exponents[1]

# Print the final explanation and result
print(f"The size of the sample space is {num_sides}^{num_dice}.")
print(f"The prime factorization of the size of the sample space is 2^{exponent_p1} * 3^{exponent_p2}.")
print("The largest possible value of m is the sum of these exponents.")
print(f"m = {exponent_p1} + {exponent_p2} = {max_m}")
