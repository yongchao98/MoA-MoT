def calculate_score(s):
    score = 0
    i = 0
    while i < len(s):
        if i < len(s) - 1 and s[i] == s[i + 1]:
            if s[i] == 'A':
                score += 1
            elif s[i] == 'B':
                score += 2
            # Skip the next character as it is part of the current sequence
            i += 1
        i += 1
    return score

# Calculate the score for the string "BBAAB"
total_score = calculate_score("BBAAB")
print(total_score)