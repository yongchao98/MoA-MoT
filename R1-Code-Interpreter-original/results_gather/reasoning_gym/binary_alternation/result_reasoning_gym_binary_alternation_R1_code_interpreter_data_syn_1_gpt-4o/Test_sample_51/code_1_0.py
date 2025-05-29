def check_feasibility(s):
    count_0 = s.count('0')
    count_1 = s.count('1')
    return count_0, count_1, abs(count_0 - count_1) > 1

# Given binary string
binary_string = "11010100111000001110100111"
count_0, count_1, is_impossible = check_feasibility(binary_string)
print(f"Count of '0's: {count_0}, Count of '1's: {count_1}, Impossible to alternate: {is_impossible}")