def get_valid_words():
    # A more comprehensive list focusing on another possible path
    return {
        'PSIS', 'PSIS', 'PITS', 'PETS', 'TETS', 'TEAT',  # New main path
        'PASS', 'PAST', 'PEST', 'PEAT', 'TEAT',          # Alternative path
        'PIES', 'TIES', 'TEAS', 'TEAT',                  # Another alternative
        'PEAS', 'TEAS', 'TEAT',                          # Yet another alternative
        'PUTS', 'POTS', 'PETS', 'TETS',                  # Additional words
        'SITS', 'SETS', 'TETS',                          # More alternatives
        'PITS', 'PETS', 'TETS', 'TEXT'                   # Extra words
    }

def one_letter_diff(word1, word2):
    return sum(c1 != c2 for c1, c2 in zip(word1, word2)) == 1

def find_all_paths(start, end, max_length=7):
    from collections import deque
    
    words = get_valid_words()
    queue = deque([(start, [start])])
    all_paths = []
    
    while queue:
        current, path = queue.popleft()
        if len(path) > max_length:
            continue
            
        if current == end:
            all_paths.append(path)
            continue
            
        for word in words:
            if word not in path and one_letter_diff(current, word):
                new_path = path + [word]
                queue.append((word, new_path))
    
    return all_paths

# Find all possible paths
paths = find_all_paths('PSIS', 'TEAT')
if paths:
    # Get the shortest path
    shortest_path = min(paths, key=len)
    print(','.join(shortest_path))
else:
    print("No solution found")