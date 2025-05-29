def one_letter_diff(word1, word2):
    diff = 0
    for c1, c2 in zip(word1, word2):
        if c1 != c2:
            diff += 1
    return diff == 1

# Potential paths to test
paths = [
    ['SHIV', 'SHIP', 'SHOP', 'SLOP', 'SLAP', 'SNAP', 'TRAP', 'TARP'],
    ['SHIV', 'SHIN', 'THIN', 'THAN', 'THAT', 'THAT', 'TRAP', 'TARP'],
    ['SHIV', 'SHIM', 'SLIM', 'SLAP', 'SNAP', 'TRAP', 'TARP']
]

def verify_path(path):
    # Check if each consecutive pair differs by exactly one letter
    for i in range(len(path)-1):
        if not one_letter_diff(path[i], path[i+1]):
            return False
    return True

# Test all paths and print the first valid one found
valid_path = None
for path in paths:
    if verify_path(path):
        valid_path = path
        print(','.join(path))
        break

if valid_path is None:
    print("No valid path found")