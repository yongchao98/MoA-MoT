from collections import deque

# More focused word list with verified words
words = '''
PLOY PLAY CLAY CLAP COMP
PLOP CLOP CHOP COMP
PLOT CLOT CLOP COMP
CAMP CHAP CHOP COMP
'''.split()

def hamming_distance(word1, word2):
    return sum(1 for a, b in zip(word1, word2) if a != b)

def find_all_paths(start, end, words):
    word_set = set(words)
    if start not in word_set or end not in word_set:
        return []
    
    paths = []
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        current = path[-1]
        
        if current == end:
            paths.append(path)
            continue
        
        for word in word_set:
            if word not in path and hamming_distance(current, word) == 1:
                new_path = path + [word]
                queue.append(new_path)
    
    return paths

paths = find_all_paths('PLOY', 'COMP', words)
for path in paths:
    print(','.join(path))