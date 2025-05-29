def calculate_score(s):
    total_score = 0
    i = 0
    while i < len(s):
        count = 1
        while i + 1 < len(s) and s[i] == s[i + 1]:
            count += 1
            i += 1
        if count >= 2:
            if s[i] == 'A':
                total_score += 1
            elif s[i] == 'B':
                total_score += 2
            elif s[i] == 'C':
                total_score += 3
            elif s[i] == 'D':
                total_score += 4
            elif s[i] == 'E':
                total_score += 5
        i += 1
    return total_score

# Given string
s = "EBAAAAACACAACAADAABC"
# Calculate and print the total score
print(calculate_score(s))