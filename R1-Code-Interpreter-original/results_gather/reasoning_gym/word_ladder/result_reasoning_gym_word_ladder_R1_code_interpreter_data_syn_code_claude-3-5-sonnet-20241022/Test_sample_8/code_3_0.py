def is_one_letter_diff(word1, word2):
    return sum(1 for a, b in zip(word1, word2) if a != b) == 1

# Define the proposed path
path = ['TOWS', 'SOWS', 'SOLS', 'SOIL']

# Verify each word differs by exactly one letter
valid = True
for i in range(len(path) - 1):
    if not is_one_letter_diff(path[i], path[i+1]):
        valid = False
        break

# Print the path if valid
if valid:
    print(','.join(path))
else:
    # Let's try an alternative path
    alt_path = ['TOWS', 'TONS', 'TOIL', 'SOIL']
    valid = True
    for i in range(len(alt_path) - 1):
        if not is_one_letter_diff(alt_path[i], alt_path[i+1]):
            valid = False
            break
    if valid:
        print(','.join(alt_path))
    else:
        print("No valid path found")