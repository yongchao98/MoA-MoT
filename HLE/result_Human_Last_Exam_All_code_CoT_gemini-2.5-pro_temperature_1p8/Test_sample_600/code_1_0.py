import sys

# Since the problem statement guarantees the last digit is the same for every positive integer k,
# we can pick any k >= 1 to find the answer. Let's use k=2 as an example.
k = 2

# According to the problem definitions:
# A_k = (10^(k+1) - 1) / 9
# B_k = 10^k
try:
    A_k = (10**(k + 1) - 1) // 9
    B_k = 10**k
except OverflowError:
    print(f"k={k} is too large to compute A_k and B_k directly. Please use a smaller k.", file=sys.stderr)
    sys.exit(1)


print(f"For k = {k}, we have:")
print(f"A_k = (10^({k}+1) - 1) / 9 = {A_k}")
print(f"B_k = 10^{k} = {B_k}")
print("-" * 30)

# We want to find the last digit of A_k^B_k - B_k^A_k.
# This is equivalent to calculating (A_k^B_k - B_k^A_k) % 10.
# The full expression is: {A_k}^{B_k} - {B_k}^{A_k}
print(f"The expression is: {A_k}^{B_k} - {B_k}^{A_k}")
print("-" * 30)

# We can find the last digit of each term separately using modular arithmetic.

# Last digit of the first term: A_k^B_k
# The numbers are too large for a direct computation, so we use the pow(base, exp, mod) function.
term1_last_digit = pow(A_k, B_k, 10)
print(f"The last digit of the first term, {A_k}^{B_k}, is computed as ({A_k} ^ {B_k}) mod 10.")
print(f"This is equivalent to ({A_k} mod 10) raised to the power of {B_k}, all modulo 10.")
print(f"The last digit of {A_k} is {A_k % 10}.")
print(f"Therefore, the last digit of {A_k}^{B_k} is {term1_last_digit}.")
print("-" * 30)


# Last digit of the second term: B_k^A_k
# A_k is guaranteed to be a positive integer for k>=1.
term2_last_digit = pow(B_k, A_k, 10)
print(f"The last digit of the second term, {B_k}^{A_k}, is computed as ({B_k} ^ {A_k}) mod 10.")
print(f"This is equivalent to ({B_k} mod 10) raised to the power of {A_k}, all modulo 10.")
print(f"The last digit of {B_k} is {B_k % 10}.")
print(f"Therefore, the last digit of {B_k}^{A_k} is {term2_last_digit}.")
print("-" * 30)


# Calculate the final result from the last digits of the two terms
# We use (a - b + 10) % 10 to handle potential negative results correctly
final_digit = (term1_last_digit - term2_last_digit + 10) % 10

print("To find the last digit of the entire expression, we subtract the last digits:")
print(f"Final last digit = (last digit of first term - last digit of second term) mod 10")
print(f"The numbers in the final equation (for the last digits) are:")
print(f"{term1_last_digit} - {term2_last_digit}  (mod 10)")
print(f"Result = {final_digit}")
print("-" * 30)

print(f"The last digit of A_k^B_k - B_k^A_k is {final_digit}.")
