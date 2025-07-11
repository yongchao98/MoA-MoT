import math

# Problem parameters
n = 99

# The sum S can be derived to be S = (1/(n+1)) * sum_{j=1 to n} [ C(n+1, j) * j^n * (-1)^(n-j) ]
# where C(n, k) is the binomial coefficient "n choose k".
# The code below calculates this sum. Python's arbitrary-precision integers handle the large numbers.

# We initialize the sum. The j=0 term is 0 since j^n = 0. We sum from j=1 to n.
sum_val = 0
for j in range(1, n + 1):
    # Calculate each component of the term
    # Binomial coefficient C(n+1, j)
    comb = math.comb(n + 1, j)
    
    # Power j^n
    j_pow_n = j**n
    
    # Sign (-1)^(n-j)
    sign = -1 if (n - j) % 2 != 0 else 1
    
    # Add the full term to our sum
    term = comb * j_pow_n * sign
    sum_val += term

# The final result is the sum divided by (n+1).
# We use integer division // as the result is known to be an exact integer.
result = sum_val // (n + 1)

# The result is (n+1)^(n-1) = 100^98 = 10^196.
# Let's present the answer as requested.
final_base = 10
final_exponent = (n - 1) * 2

print("The calculation is based on a simplified summation formula derived from the problem statement.")
print(f"For n = {n}, the result is ({n}+1)**({n}-1).")
print(f"This is equal to {n+1}**{n-1}, which can be written as a power of 10.")
# Fulfilling the request: "output each number in the final equation!"
print(f"{final_base} ** {final_exponent} = {result}")
