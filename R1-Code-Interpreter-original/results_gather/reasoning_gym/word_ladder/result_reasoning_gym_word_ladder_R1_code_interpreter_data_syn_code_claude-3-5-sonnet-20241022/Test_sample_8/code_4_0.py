def is_one_letter_diff(word1, word2):
    diff_count = 0
    for c1, c2 in zip(word1, word2):
        if c1 != c2:
            diff_count += 1
    return diff_count == 1

# Let's try this path
path = ['TOWS', 'TOWN', 'TORN', 'TOIL', 'SOIL']

# Check if each consecutive pair differs by exactly one letter
# and print the differences for analysis
for i in range(len(path)-1):
    w1, w2 = path[i], path[i+1]
    print(f"Checking {w1} -> {w2}:")
    diff_positions = [(j, w1[j], w2[j]) for j in range(4) if w1[j] != w2[j]]
    print(f"Differences: {diff_positions}")
    print(f"Valid step: {is_one_letter_diff(w1, w2)}")

# Verify the entire path
valid = all(is_one_letter_diff(path[i], path[i+1]) for i in range(len(path)-1))
print("\nFinal path:", ','.join(path) if valid else "No valid path found")