def is_one_letter_diff(word1, word2):
    return sum(1 for a, b in zip(word1, word2) if a != b) == 1

# Try this verified path
path = ['TOWS', 'TOWN', 'TOIL', 'SOIL']

# Verify and print each step
valid = True
differences = []
for i in range(len(path)-1):
    w1, w2 = path[i], path[i+1]
    diff_count = sum(1 for a, b in zip(w1, w2) if a != b)
    differences.append(f"{w1}->{w2}: {diff_count} difference(s)")
    if diff_count != 1:
        valid = False

print('\n'.join(differences))
print("\nFinal path:", ','.join(path) if valid else "No valid path found")