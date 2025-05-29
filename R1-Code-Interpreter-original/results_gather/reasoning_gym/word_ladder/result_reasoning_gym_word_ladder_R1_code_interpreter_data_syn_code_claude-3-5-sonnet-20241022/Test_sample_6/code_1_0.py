from collections import deque

# Expanded list of common 4-letter words
words = '''
PLOY PLAY CLAY CLAW CLAP CHAP CHOP COMP
PLOT PLAN CLAN CLAM CRAM CRAP CROP
PLOD CLOD CLOT COAT CHAT CHIP
PLOW CLOW CLOG COMB
PLUG CLUE COPE
PLUS CLUE COME
PLOP CLOP
CAMP DAMP DAME DOME COME COMP
CHAP CHIP CHOP COMP
PLOY PLOT CLOT CLOP COMP
PLOY PLAY CLAY CLAP CHAP COMP
'''.split()

def get_neighbors(word, word_set):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in word_set:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end, words):
    word_set = set(words)
    if start not in word_set or end not in word_set:
        return None
    
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        current = path[-1]
        
        if current == end:
            return path
        
        for next_word in get_neighbors(current, word_set):
            if next_word not in visited:
                visited.add(next_word)
                queue.append(path + [next_word])
    
    return None

# Try multiple potential paths
paths = [
    find_path('PLOY', 'COMP', words),
]

# Print all found paths
for path in paths:
    if path:
        print(','.join(path))