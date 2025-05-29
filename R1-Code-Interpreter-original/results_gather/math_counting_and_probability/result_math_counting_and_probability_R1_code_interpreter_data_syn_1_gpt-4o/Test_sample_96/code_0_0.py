# Calculate the number of three-digit numbers that are multiples of 5
multiples_of_5 = (995 - 100) // 5 + 1

# Calculate the number of three-digit numbers that are multiples of 7
multiples_of_7 = (994 - 105) // 7 + 1

# Calculate the number of three-digit numbers that are multiples of both 5 and 7 (i.e., multiples of 35)
multiples_of_35 = (980 - 105) // 35 + 1

# Apply the principle of inclusion-exclusion
multiples_of_5_or_7 = multiples_of_5 + multiples_of_7 - multiples_of_35

# Total number of three-digit numbers
total_three_digit_numbers = 999 - 100 + 1

# Calculate the number of three-digit numbers that are multiples of neither 5 nor 7
multiples_of_neither = total_three_digit_numbers - multiples_of_5_or_7

print(multiples_of_neither)