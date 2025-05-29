def get_valid_words():
    # A carefully curated set of words specifically for this puzzle
    return {
        'PSIS', 'PASS', 'PAST', 'PEST', 'PEAT', 'TEAT',  # Main path
        'TEAS', 'PEAS', 'PSST', 'POST', 'TEST',          # Alternative words
    }

def one_letter_diff(word1, word2):
    # Check if exactly one letter is different between two words
    diff = 0
    for c1, c2 in zip(word1, word2):
        if c1 != c2:
            diff += 1
        if diff > 1:
            return False
    return diff == 1

def find_all_paths(current, end, path, all_paths, words, used):
    if current == end:
        all_paths.append(path[:])
        return
    
    for word in words:
        if word not in used and one_letter_diff(current, word):
            used.add(word)
            path.append(word)
            find_all_paths(word, end, path, all_paths, words, used)
            path.pop()
            used.remove(word)

# Find all possible paths
start = 'PSIS'
end = 'TEAT'
words = get_valid_words()
all_paths = []
used = {start}
find_all_paths(start, end, [start], all_paths, words, used)

# Print the shortest valid path if found
if all_paths:
    shortest_path = min(all_paths, key=len)
    print(','.join(shortest_path))
else:
    print("No solution found")