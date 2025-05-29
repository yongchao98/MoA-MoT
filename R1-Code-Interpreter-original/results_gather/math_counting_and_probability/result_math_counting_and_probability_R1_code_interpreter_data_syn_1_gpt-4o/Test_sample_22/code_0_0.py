# Total number of three-digit numbers
total_numbers = 999 - 100 + 1

# Number of three-digit numbers without any 7 or 9
without_7_or_9 = 7 * 8 * 8

# Number of three-digit numbers with at least one 7 or 9
with_7_or_9 = total_numbers - without_7_or_9

print(with_7_or_9)