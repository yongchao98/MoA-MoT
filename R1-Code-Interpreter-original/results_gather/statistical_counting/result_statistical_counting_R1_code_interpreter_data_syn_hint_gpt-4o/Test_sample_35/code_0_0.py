def calculate_score(s):
    total_score = 0
    i = 0
    while i < len(s) - 1:
        if s[i] == s[i + 1]:
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
            elif s[i] == 'F':
                total_score += 6
            # Move the index to skip the consecutive characters
            while i < len(s) - 1 and s[i] == s[i + 1]:
                i += 1
        i += 1
    return total_score

# Given string
s = "FADAADABCEBCDFBAFDFAAFFFE"
# Calculate and print the total score
print(calculate_score(s))