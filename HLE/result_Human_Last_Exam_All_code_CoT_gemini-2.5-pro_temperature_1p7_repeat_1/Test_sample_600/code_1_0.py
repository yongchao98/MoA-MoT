# The user wants to find the last digit of A_k^(B_k) - B_k^(A_k) for any positive integer k.
# The problem states that this last digit is always the same.
# We can therefore pick any positive integer for k, for example k=1, to find the answer.

k = 1

# Calculate A_k and B_k for the chosen k.
# A_k is a repunit of length k+1.
# Note: Using integer division // to ensure the result is an integer.
A_k = (10**(k + 1) - 1) // 9
# B_k is 10 to the power of k.
B_k = 10**k

# To find the last digit of the expression, we use modular arithmetic (modulo 10).

# Find the last digit of the first term: A_k ^ B_k
# The last digit of A_k is A_k % 10.
last_digit_Ak = A_k % 10
# Any integer power of a number ending in 1 will also end in 1.
# B_k is a positive integer for k>=1.
ld_term1 = pow(last_digit_Ak, B_k, 10)

# Find the last digit of the second term: B_k ^ A_k
# The last digit of B_k is B_k % 10.
last_digit_Bk = B_k % 10
# Any positive integer power of a number ending in 0 is 0.
# A_k is a positive integer for k>=1.
ld_term2 = pow(last_digit_Bk, A_k, 10)

# The last digit of the difference is (ld_term1 - ld_term2) % 10.
# We add 10 before the modulo to ensure the result is non-negative.
final_last_digit = (ld_term1 - ld_term2 + 10) % 10

# Print the analysis for the chosen k as requested.
print(f"For k = {k}:")
print(f"A_k = (10^({k}+1) - 1) / 9 = {A_k}")
print(f"B_k = 10^{k} = {B_k}")
print(f"The equation is {A_k}^{B_k} - {B_k}^{A_k}.")
print(f"The last digit of the first term, {A_k}^{B_k}, is {ld_term1}.")
print(f"The last digit of the second term, {B_k}^{A_k}, is {ld_term2}.")
print(f"The last digit of the expression {A_k}^{B_k} - {B_k}^{A_k} is the last digit of a number ending in {ld_term1} minus a number ending in {ld_term2}.")
print(f"This is calculated as ({ld_term1} - {ld_term2}) mod 10, which is {final_last_digit}.")
print(f"\nThe constant last digit is {final_last_digit}.")
