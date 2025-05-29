def check_digit(ch1, ch2):
    row1, col1 = (ch1 - 1) // 3, (ch1 - 1) % 3
    row2, col2 = (ch2 - 1) // 3, (ch2 - 1) % 3
    return abs(row1 - row2) + abs(col1 - col2)

def security_check(security_code, user_input):
    if len(security_code) != len(user_input):
        return False
    count_error = 0
    for i in range(len(security_code)):
        if security_code[i] != user_input[i]:
            if check_digit(int(security_code[i]), int(user_input[i])) != 1:
                return False
            count_error += 1
    return count_error <= 1

# Given input
security_code = '923745'
user_input = '623745'

# Perform the security check
result = security_check(security_code, user_input)

# Print the result
print(result)