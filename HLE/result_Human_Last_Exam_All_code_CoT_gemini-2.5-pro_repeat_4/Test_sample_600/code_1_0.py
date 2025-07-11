import sys

# The problem states the last digit is the same for every positive integer k.
# We will demonstrate the calculation for a sample value, k=2.
k = 2

# Define A_k and B_k.
# A_k is a number consisting of k+1 ones.
A_k = (10**(k+1) - 1) // 9
# B_k is 10 to the power of k.
B_k = 10**k

# The task is to find the last digit of A_k**B_k - B_k**A_k.
# This is equivalent to calculating (A_k**B_k - B_k**A_k) % 10.
# We can use modular arithmetic to simplify the calculation.

print(f"This script calculates the last digit of A_k^B_k - B_k^A_k for a given k.")
print(f"The result is constant for all positive integers k. We will use k={k} as an example.")
print("-" * 30)

# Print the definitions and values for the chosen k
print(f"A_k is defined as (10^(k+1) - 1) / 9.")
print(f"For k={k}, A_{k} = (10^({k}+1) - 1) / 9 = {A_k}.")
print(f"B_k is defined as 10^k.")
print(f"For k={k}, B_{k} = 10^{k} = {B_k}.")

# Print the specific expression for k
print("\nThe expression for k=2 is:")
print(f"{A_k}^{B_k} - {B_k}^{A_k}")
print("-" * 30)

# Step 1: Find the last digit of the first term, A_k**B_k.
# The last digit of a number 'n' is 'n % 10'.
# The last digit of n^p is pow(n, p, 10).
last_digit_A_k = A_k % 10
# The last digit of any A_k is 1.
# Any power of a number ending in 1 also ends in 1.
term1_last_digit = pow(last_digit_A_k, B_k, 10)

print("Step 1: Find the last digit of the first term.")
print(f"The last digit of A_k ({A_k}) is {last_digit_A_k}.")
print(f"The last digit of A_k^B_k is the same as the last digit of {last_digit_A_k}^{B_k}.")
print(f"The last digit of the first term is {term1_last_digit}.")
print("-" * 30)

# Step 2: Find the last digit of the second term, B_k**A_k.
last_digit_B_k = B_k % 10
# The last digit of any B_k (for k>=1) is 0.
# Any positive power of a number ending in 0 also ends in 0.
# A_k is always a positive integer for k>=1.
if A_k == 0:
    # This case does not happen for k>=1, but it's good practice for pow(0,0)
    term2_last_digit = 1
else:
    term2_last_digit = pow(last_digit_B_k, A_k, 10)

print("Step 2: Find the last digit of the second term.")
print(f"The last digit of B_k ({B_k}) is {last_digit_B_k}.")
print(f"The last digit of B_k^A_k is the same as the last digit of {last_digit_B_k}^{A_k}.")
print(f"The last digit of the second term is {term2_last_digit}.")
print("-" * 30)

# Step 3: Find the last digit of the difference.
# The result is (term1_last_digit - term2_last_digit) % 10.
# We add 10 before the modulo to ensure the result is non-negative.
final_last_digit = (term1_last_digit - term2_last_digit + 10) % 10

print("Step 3: Combine the results.")
print("The last digit of the expression is the last digit of the difference of the two terms' last digits.")
print(f"Last digit of ({term1_last_digit} - {term2_last_digit}) = ({term1_last_digit - term2_last_digit})")
print(f"({term1_last_digit - term2_last_digit}) mod 10 = {final_last_digit}")
print("-" * 30)

print(f"The last digit of A_k^B_k - B_k^A_k for any positive integer k is {final_last_digit}.")
