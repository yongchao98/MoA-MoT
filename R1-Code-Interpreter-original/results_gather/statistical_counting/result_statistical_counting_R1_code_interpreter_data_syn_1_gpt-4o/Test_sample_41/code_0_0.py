def calculate_score(s):
    total_score = 0
    current_char = ''
    count = 0

    for char in s:
        if char == current_char:
            count += 1
        else:
            if count >= 2:
                if current_char == 'A':
                    total_score += 1
                elif current_char == 'B':
                    total_score += 2
                elif current_char == 'C':
                    total_score += 3
            current_char = char
            count = 1

    # Check the last sequence
    if count >= 2:
        if current_char == 'A':
            total_score += 1
        elif current_char == 'B':
            total_score += 2
        elif current_char == 'C':
            total_score += 3

    return total_score

# Given string
s = "AAAACCCCBCAAABCACCBAABAAAABACBCABAA"
# Calculate and print the total score
print(calculate_score(s))