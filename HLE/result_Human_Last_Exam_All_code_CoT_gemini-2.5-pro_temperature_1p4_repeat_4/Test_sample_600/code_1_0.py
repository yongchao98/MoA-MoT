# We are asked to find the last digit of A_k^(B_k) - B_k^(A_k)
# where A_k = (10**(k+1)-1)/9 and B_k = 10**k.
# The problem states the last digit is the same for every positive integer k.
# We can choose any positive integer for k to find the answer, let's use k=2.

k = 2

# Calculate A_k and B_k for k=2
# A_k will be a number made of k+1 ones. For k=2, A_k is 111.
A_k = (10**(k + 1) - 1) // 9
# B_k is 10 to the power of k. For k=2, B_k is 100.
B_k = 10**k

print(f"For k={k}, we are finding the last digit of the equation:")
print(f"{A_k}^{B_k} - {B_k}^{A_k}")
print("-" * 30)

# To find the last digit, we use the modulo 10 operator.
# The last digit of X^Y can be found efficiently with pow(X, Y, 10).

# Find the last digit of the first term: A_k^B_k
# The last digit of A_k (111) is 1. Any integer power of a number ending in 1 also ends in 1.
last_digit_term1 = pow(A_k, B_k, 10)

# Find the last digit of the second term: B_k^A_k
# The last digit of B_k (100) is 0. Any positive integer power of a number ending in 0 also ends in 0.
# A_k (111) is a positive integer.
last_digit_term2 = pow(B_k, A_k, 10)

# The last digit of the expression is (last_digit_term1 - last_digit_term2) % 10.
# We add 10 before the modulo to handle potential negative results correctly (e.g., (3 - 5) % 10).
final_last_digit = (last_digit_term1 - last_digit_term2 + 10) % 10

print(f"The last digit of {A_k}^{B_k} is {last_digit_term1}.")
print(f"The last digit of {B_k}^{A_k} is {last_digit_term2}.")
print(f"The last digit of the expression is ({last_digit_term1} - {last_digit_term2}) mod 10, which is {final_last_digit}.")
