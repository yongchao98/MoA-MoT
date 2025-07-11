# We want to find the last digit of A_k^(B_k) - B_k^(A_k) for any positive integer k.
# This is equivalent to calculating (A_k^(B_k) - B_k^(A_k)) mod 10.

# For any positive integer k, A_k is a number of the form 11...1.
# Its last digit is always 1.
last_digit_Ak = 1

# For any positive integer k, B_k is a power of 10 (10, 100, ...).
# Its last digit is always 0.
last_digit_Bk = 0

# Now we find the last digit of A_k^(B_k).
# This is determined by the last digit of A_k.
# For k>=1, B_k is a positive integer.
# The last digit of a number ending in 1 raised to any positive integer power is 1.
last_digit_Ak_pow_Bk = 1

# Next, we find the last digit of B_k^(A_k).
# This is determined by the last digit of B_k.
# For k>=1, A_k is a positive integer.
# The last digit of a number ending in 0 raised to any positive integer power is 0.
last_digit_Bk_pow_Ak = 0

# The last digit of the difference is the result of the subtraction of the last digits.
# We use modulo 10 to handle the result properly.
# (a - b) mod n = (a mod n - b mod n + n) mod n
final_last_digit = (last_digit_Ak_pow_Bk - last_digit_Bk_pow_Ak + 10) % 10

print("The problem is to find the last digit of A_k^B_k - B_k^A_k.")
print("The last digit of A_k^B_k is the last digit of (...1) raised to a positive power, which is 1.")
print("The last digit of B_k^A_k is the last digit of (...0) raised to a positive power, which is 0.")
print("The final equation for the last digit is:")
print(f"{last_digit_Ak_pow_Bk} - {last_digit_Bk_pow_Ak} = {final_last_digit}")
print(f"Thus, the last digit is always {final_last_digit}.")
