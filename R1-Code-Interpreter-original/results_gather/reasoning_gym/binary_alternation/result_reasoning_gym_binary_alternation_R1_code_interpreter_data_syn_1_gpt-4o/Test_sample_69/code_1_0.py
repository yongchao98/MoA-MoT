def check_counts(s):
    count_0 = s.count('0')
    count_1 = s.count('1')
    return count_0, count_1, abs(count_0 - count_1)

# Given binary string
binary_string = "010001011100011101001101001010"
counts = check_counts(binary_string)
print(counts)