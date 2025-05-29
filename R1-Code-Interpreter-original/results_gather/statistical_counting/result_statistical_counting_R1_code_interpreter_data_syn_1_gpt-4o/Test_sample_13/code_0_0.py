def calculate_score(s):
    score = 0
    i = 0
    while i < len(s) - 1:
        if s[i] == s[i + 1]:
            if s[i] == 'A':
                score += 1
            elif s[i] == 'B':
                score += 2
            elif s[i] == 'C':
                score += 3
            elif s[i] == 'D':
                score += 4
            elif s[i] == 'E':
                score += 5
            # Move to the next different character
            while i < len(s) - 1 and s[i] == s[i + 1]:
                i += 1
        i += 1
    return score

# Given string
s = "CBABBECAECEEBBCEECBB"
# Calculate and print the total score
print(calculate_score(s))