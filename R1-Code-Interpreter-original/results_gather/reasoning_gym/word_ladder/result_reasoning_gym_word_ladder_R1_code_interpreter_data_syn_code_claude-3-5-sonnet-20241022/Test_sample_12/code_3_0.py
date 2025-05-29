def get_valid_words():
    # Absolutely minimal set of words needed for the path
    return {
        'SNOG',
        'SNUG',
        'SMUG',
        'SLUG',
        'SLOG',
        'SAGS',
        'SOGS',
        'SUGS',
        'SUQS'
    }

def is_one_letter_different(word1, word2):
    differences = 0
    for c1, c2 in zip(word1, word2):
        if c1 != c2:
            differences += 1
    return differences == 1

def find_all_paths(current, end, path, valid_words, all_paths):
    if current == end:
        all_paths.append(path)
        return
    
    for word in valid_words:
        if word not in path and is_one_letter_different(current, word):
            find_all_paths(word, end, path + [word], valid_words, all_paths)

all_paths = []
valid_words = get_valid_words()
find_all_paths('SNOG', 'SUQS', ['SNOG'], valid_words, all_paths)

if all_paths:
    # Print the shortest path found
    shortest_path = min(all_paths, key=len)
    print(','.join(shortest_path))
else:
    print("No valid path found")