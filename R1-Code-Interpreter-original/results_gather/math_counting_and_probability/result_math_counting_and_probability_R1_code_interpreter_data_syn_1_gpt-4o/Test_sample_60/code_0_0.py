# Total number of three-digit numbers without a '5'
total_without_5 = 8 * 9 * 9

# Total number of three-digit numbers without a '5' and without a '3'
total_without_5_and_3 = 7 * 8 * 8

# Total number of three-digit numbers with at least one '3' and no '5'
result = total_without_5 - total_without_5_and_3

print(result)