# The problem states that the last digit of the expression is the same for every positive integer k.
# We can therefore choose any value for k to find the result. We will use k=1 as a demonstration.

k = 1

# Step 1: Define A_k and B_k for our chosen k.
# A_k = (10**(k+1) - 1) / 9. For k=1, A_1 = (10^2 - 1)/9 = 99/9 = 11.
# B_k = 10**k. For k=1, B_1 = 10^1 = 10.
A_k_val = (10**(k + 1) - 1) // 9
B_k_val = 10**k

# Step 2: Find the last digit of the first term, A_k^(B_k).
# The last digit of a number n is n % 10.
# The last digit of A_k is always 1 (since A_k is 11, 111, 1111, ...).
# Any positive integer power of a number ending in 1 also ends in 1.
# We can confirm this using modular exponentiation: pow(base, exp, mod).
last_digit_term1 = pow(A_k_val, B_k_val, 10)

# Step 3: Find the last digit of the second term, B_k^(A_k).
# The last digit of B_k (10^k) is always 0 for k>=1.
# Any positive integer power of a number ending in 0 (like 10, 20, etc.) also ends in 0.
# pow(0, exp, mod) is 0 for exp > 0.
last_digit_term2 = pow(B_k_val, A_k_val, 10)

# Step 4: Find the last digit of the difference A_k^(B_k) - B_k^(A_k).
# This is equivalent to (last_digit_term1 - last_digit_term2) mod 10.
# We add 10 before taking the modulo to ensure the result is positive.
final_last_digit = (last_digit_term1 - last_digit_term2 + 10) % 10

# Step 5: Print the results, showing each number in the final equation.
print(f"For k = {k}:")
print(f"A_{k} = {A_k_val}")
print(f"B_{k} = {B_k_val}")
print(f"We want the last digit of the expression: {A_k_val}^{B_k_val} - {B_k_val}^{A_k_val}")
print(f"The last digit of the first term, {A_k_val}^{B_k_val}, is {last_digit_term1}.")
print(f"The last digit of the second term, {B_k_val}^{A_k_val}, is {last_digit_term2}.")
print(f"The last digit of the difference corresponds to ({last_digit_term1} - {last_digit_term2}) mod 10.")
print(f"Final calculated last digit: {final_last_digit}")