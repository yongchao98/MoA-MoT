import math

# The problem states the last digit is the same for every positive integer k.
# We will use k=1 for this demonstration.
k = 1

# Calculate A_k for the chosen k.
# A_k is a sequence of (k+1) ones.
A_k = (10**(k+1) - 1) // 9

# Calculate B_k for the chosen k.
# B_k is 10 to the power of k.
B_k = 10**k

print(f"This script calculates the last digit of the expression A_k^(B_k) - B_k^(A_k).")
print(f"The result is the same for any positive integer k. We'll use k = {k} as an example.")
print("-" * 30)
print(f"For k = {k}:")
print(f"A_{k} = (10^({k+1}) - 1) / 9 = {A_k}")
print(f"B_{k} = 10^{k} = {B_k}")

print(f"\nThe expression is: {A_k}^{B_k} - {B_k}^{A_k}")
print("\nTo find the last digit, we compute the expression modulo 10.")

# To find the last digit of A_k^(B_k), we can compute (A_k % 10) ^ B_k.
# Since A_k always ends in 1, any power of it will also end in 1.
term1_last_digit = pow(A_k, B_k, 10)
print(f"The last digit of the first term ({A_k}^{B_k}) is {term1_last_digit}.")

# To find the last digit of B_k^(A_k), we can compute (B_k % 10) ^ A_k.
# Since B_k is a multiple of 10, it ends in 0. Any positive power of it will also end in 0.
term2_last_digit = pow(B_k, A_k, 10)
print(f"The last digit of the second term ({B_k}^{A_k}) is {term2_last_digit}.")

# The final calculation for the last digit is the last digit of the difference.
# The addition of 10 before the modulo handles potential negative results (e.g., 3 - 7).
final_last_digit = (term1_last_digit - term2_last_digit + 10) % 10

print(f"\nThe final equation for the last digit is based on the last digits of the terms:")
print(f"({term1_last_digit} - {term2_last_digit}) mod 10")

print(f"\nTherefore, the last digit of the expression is {final_last_digit}.")