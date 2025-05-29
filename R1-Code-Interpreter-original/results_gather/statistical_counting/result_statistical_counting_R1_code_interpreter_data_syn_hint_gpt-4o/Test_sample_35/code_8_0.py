def calculate_score(s):
    score_map = {'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5, 'F': 6}
    total_score = 0
    i = 0
    while i < len(s) - 1:
        if s[i] == s[i + 1]:
            total_score += score_map[s[i]]
            # Move the index to skip the consecutive characters
            while i < len(s) - 1 and s[i] == s[i + 1]:
                i += 1
        i += 1
    return total_score

# Given string
s = "FADAADABCEBCDFBAFDFAAFFFE"
# Calculate and print the total score
print(calculate_score(s))