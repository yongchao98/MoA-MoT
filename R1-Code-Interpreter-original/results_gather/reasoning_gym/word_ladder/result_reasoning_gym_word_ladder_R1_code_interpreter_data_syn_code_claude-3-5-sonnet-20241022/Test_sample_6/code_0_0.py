from collections import deque, defaultdict

# Common 4-letter words
words = '''PLOY PLAY CLAY CLAW CLAP CHAP CHOP COMP
PLOT PLAN CLAN CLAM CRAM CRAP CROP
PLOD CLOD CLOT COAT CHAT CHIP
PLOW CLOW CLOG COMB
PLUG CLUE COPE
PLUS CLUE COME'''.split()

def get_neighbors(word):
    neighbors = []
    for i in range(4):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in word_set and new_word != word:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    if start not in word_set or end not in word_set:
        return None
    
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        word = path[-1]
        
        if word == end:
            return path
            
        for neighbor in get_neighbors(word):
            if neighbor not in visited:
                visited.add(neighbor)
                new_path = list(path)
                new_path.append(neighbor)
                queue.append(new_path)
    
    return None

word_set = set(words)
path = find_path('PLOY', 'COMP')
if path:
    print(','.join(path))
else:
    print("No path found")