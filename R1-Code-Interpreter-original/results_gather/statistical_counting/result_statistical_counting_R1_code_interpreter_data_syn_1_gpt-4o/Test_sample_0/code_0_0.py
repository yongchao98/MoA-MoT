def calculate_score(s):
    total_score = 0
    i = 0
    while i < len(s) - 1:
        if s[i] == s[i + 1]:
            if s[i] == 'A':
                total_score += 1
            elif s[i] == 'B':
                total_score += 2
            # Move the index to the end of the current sequence
            while i < len(s) - 1 and s[i] == s[i + 1]:
                i += 1
        i += 1
    return total_score

# Calculate the score for the string "BBAAB"
score = calculate_score("BBAAB")
print(score)