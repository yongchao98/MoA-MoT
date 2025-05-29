def get_valid_words():
    # A different set of words focusing on another possible path
    return {
        'PSIS', 'PIES', 'TIES', 'TEAS', 'TEAT',  # New main path
        'PASS', 'PAST', 'PEST', 'PEAT',          # Alternative words
        'TESS', 'TEST', 'TEAM', 'TEAR',          # T-words
        'PEAS', 'PETS', 'PENS', 'TENS'           # P-words
    }

def one_letter_diff(word1, word2):
    diff = sum(1 for c1, c2 in zip(word1, word2) if c1 != c2)
    return diff == 1

def find_path(start, end):
    from collections import deque
    
    words = get_valid_words()
    queue = deque([(start, [start])])
    seen = {start}
    
    while queue:
        current, path = queue.popleft()
        if current == end:
            return path
        
        for word in words:
            if word not in seen and one_letter_diff(current, word):
                seen.add(word)
                new_path = path + [word]
                queue.append((word, new_path))
    return None

# Try to find the path
path = find_path('PSIS', 'TEAT')
if path:
    print(','.join(path))
else:
    print("No solution found")