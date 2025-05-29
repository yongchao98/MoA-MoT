def one_letter_diff(word1, word2):
    diff = 0
    for c1, c2 in zip(word1, word2):
        if c1 != c2:
            diff += 1
    return diff == 1

# Test a single specific path
path = ['SHIV', 'SHOW', 'SHOP', 'STOP', 'STEP', 'TRAP', 'TARP']

def verify_and_explain_path(path):
    # Check if each consecutive pair differs by exactly one letter
    # and explain the changes
    for i in range(len(path)-1):
        word1, word2 = path[i], path[i+1]
        if not one_letter_diff(word1, word2):
            print(f"Invalid transition between {word1} and {word2}")
            return False
        else:
            # Find which letter changed
            for j in range(4):
                if word1[j] != word2[j]:
                    print(f"Changed {word1} to {word2} by replacing {word1[j]} with {word2[j]}")
    return True

if verify_and_explain_path(path):
    print(','.join(path))
else:
    print("No valid path found")