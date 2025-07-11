import math

def is_prime(n):
    """
    A simple function to check if a number is prime.
    """
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_factors(n):
    """
    Finds all pairs of factors (a, b) of n such that a > b > 0.
    """
    factors = []
    # Iterate from 1 up to the square root of n
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            # i is a factor, find the other factor j
            j = n // i
            # We need a > b, so we set a = j and b = i
            if j > i:
                factors.append((j, i))
    return factors

# The number from the equation x^2 - y^2 = N
N = 2023

print(f"We are solving the equation x^2 - y^2 = {N} for positive, coprime integers (x, y).")
print(f"The equation can be factored as (x + y)(x - y) = {N}.")
print("Let a = x + y and b = x - y. We need to find factors (a, b) of 2023 such that a > b > 0.")
print("The additional conditions are that x and y must be coprime, and both a and b must be prime numbers.\n")

# Step 1: Factor 2023 into pairs (a, b)
factor_pairs = find_factors(N)
print(f"The factor pairs (a, b) of {N} with a > b are: {factor_pairs}\n")

valid_pair_count = 0
found_pairs = []

# Step 2 & 3: Iterate through pairs, compute x and y, and verify all conditions
for a, b in factor_pairs:
    print(f"--- Checking Case for (a, b) = ({a}, {b}) ---")
    
    # Verify prime condition for a and b
    a_is_prime = is_prime(a)
    b_is_prime = is_prime(b)
    print(f"Step A: Verify that a and b are prime numbers.")
    print(f"Is a = {a} prime? {a_is_prime}")
    print(f"Is b = {b} prime? {b_is_prime}")
    
    # If this condition fails, we can stop checking this pair.
    if not (a_is_prime and b_is_prime):
        print("Verification failed: Both a and b must be prime. This case is invalid.\n")
        continue

    # This part of the code will only run if a and b are both prime.
    # Step B: Compute x and y
    x = (a + b) // 2
    y = (a - b) // 2
    print(f"Step B: Compute x and y from a and b.")
    print(f"x = (a + b) / 2 = ({a} + {b}) / 2 = {x}")
    print(f"y = (a - b) / 2 = ({a} - {b}) / 2 = {y}")

    # Verify coprime condition for x and y
    gcd_xy = math.gcd(x, y)
    are_coprime = (gcd_xy == 1)
    print(f"Step C: Verify that x and y are coprime.")
    print(f"The greatest common divisor of x={x} and y={y} is {gcd_xy}.")
    print(f"Are x and y coprime? {are_coprime}")

    if are_coprime:
        print(f"Verification successful: The pair (x, y) = ({x}, {y}) satisfies all conditions.\n")
        valid_pair_count += 1
        found_pairs.append((x,y))
    else:
        print(f"Verification failed: x and y are not coprime. This case is invalid.\n")

print("--- Final Result ---")
print(f"Based on the analysis, there are {valid_pair_count} pairs (x,y) that satisfy all given conditions.")
if valid_pair_count > 0:
    for x, y in found_pairs:
        # Display the final equation for each valid pair found
        print(f"Valid pair: (x, y) = ({x}, {y}). Equation: {x}^2 - {y}^2 = {x*x} - {y*y} = {x*x - y*y}.")

print(f"<<<{valid_pair_count}>>>")