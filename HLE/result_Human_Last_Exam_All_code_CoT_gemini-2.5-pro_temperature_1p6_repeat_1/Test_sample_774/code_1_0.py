import math

def sieve_for_primes_to_n(n):
    """
    Generates primes up to n using the Sieve of Eratosthenes.
    Returns a list of primes.
    """
    size = n // 2
    sieve = [1] * size
    limit = int(n**0.5)
    for i in range(1, (limit + 1) // 2):
        if sieve[i]:
            val = 2 * i + 1
            tmp = val**2
            start = tmp // 2
            sieve[start::val] = [0] * ((size - 1 - start) // val + 1)
    return [2] + [2 * i + 1 for i, v in enumerate(sieve) if v and i > 0]

# Step 1: Analyze the sum term.
# The injectivity radius l(n,p) of the Stiefel manifold M(n,p) is known to be pi.
# The sum has 10*10 = 100 terms, so its value is 100 * pi.
sum_val_expr = "100 * pi"
sum_val_numeric = 100 * math.pi

# Step 2: Analyze the integral term.
# The integral can be simplified to I = I_1 + I_2, where I_2 is the integral of x*e^(-x) dx, which is 1.
# I_1 involves manifold dimensions. Let's calculate them.

# We need primes up to index 10231.
# An estimate for p_n is n*log(n). For n=10231, p_n ~ 10231*ln(10231) ~ 94,400.
# A sieve limit of 110,000 is safe.
sieve_limit = 110000
primes_list = sieve_for_primes_to_n(sieve_limit)

p_781 = primes_list[781 - 1]
p_2321 = primes_list[2321 - 1]
p_8231 = primes_list[8231 - 1]
p_10231 = primes_list[10231 - 1]

# Calculate dimensions d1 and d2
n1, p1 = p_8231, p_781
dim1 = n1 * p1 - (p1 * (p1 + 1)) / 2

n2, p2 = p_10231, p_2321
dim2 = n2 * p2 - (p2 * (p2 + 1)) / 2

# Based on mathematical analysis (Dominated Convergence Theorem), the integral I_1 is 0
# because the dimensions d1 and d2 are very large.
integral_val = 1.0

# Step 3: Final Calculation.
final_result = sum_val_numeric * integral_val

print("This problem evaluates a mathematical expression involving geometry and calculus.")
print("The final result is a product of two terms.")
print("\n--- Term 1: The Sum ---")
print(f"The sum is composed of 100 terms, each evaluating to pi.")
print(f"Value of the sum = {sum_val_expr} ≈ {sum_val_numeric:.4f}")

print("\n--- Term 2: The Integral ---")
print("The integral splits into two parts: I = I_1 + I_2.")
print("The second part, I_2 = integral(x*exp(-x)), evaluates to 1.")
print("The first part, I_1, depends on manifold dimensions d1 and d2.")
print("\nCalculating the prime numbers:")
print(f"p_781 = {p_781}")
print(f"p_8231 = {p_8231}")
print(f"p_2321 = {p_2321}")
print(f"p_10231 = {p_10231}")
print("\nCalculating the dimensions:")
print(f"d1 = dim(M({n1}, {p1})) = {n1}*{p1} - {p1}*({p1}+1)/2 = {int(dim1)}")
print(f"d2 = dim(M({n2}, {p2})) = {n2}*{p2} - {p2}*({p2}+1)/2 = {int(dim2)}")
print("\nBecause d1 and d2 are very large, the integral I_1 evaluates to 0.")
print(f"Value of the integral = I_1 + I_2 = 0 + 1 = {integral_val}")

print("\n--- Final Calculation ---")
print(f"Final Expression = (Sum) * (Integral)")
print(f"Result = ({sum_val_expr}) * ({integral_val})")
print(f"Final numeric value: {final_result}")