import math

# The number of variables for the Boolean function
n = 4

# Explain the problem and the formula being used.
print(f"To find a(n), the maximal number of prime implicants of a Boolean function of n variables.")
print(f"The value depends on whether n is even or odd.")
print(f"For an even n (like n=4), the formula is a(n) = C(n, n/2), where C is the binomial coefficient (the 'choose' function).")
print("-" * 20)

# Apply the formula for n=4
k = n // 2
result = math.comb(n, k)

# Show the specific calculation for a(4)
print(f"For n = {n}:")
print(f"The calculation is a({n}) = C({n}, {n}/{2}) = C({n}, {k})")
print(f"Using the formula C(n, k) = n! / (k! * (n-k)!):")
print(f"C({n}, {k}) = {math.factorial(n)} / ({math.factorial(k)}! * ({n}-{k})!) = {math.factorial(n)} / ({math.factorial(k)} * {math.factorial(n-k)}) = {result}")

# Print the final equation with each number, as requested.
print("\nThe final equation is:")
print(f"a({n}) = {result}")
