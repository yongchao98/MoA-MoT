from collections import deque, defaultdict

# First, let's define a set of common 4-letter words
words = set([
    'MAGE', 'MADE', 'MAKE', 'MALE', 'MARE', 'PAGE', 'PALE', 'PANE', 'PANG',
    'PARE', 'PARK', 'PATE', 'PAVE', 'WAGE', 'WAKE', 'WANE', 'WANG', 'WARE',
    'WARM', 'WART', 'WAVE', 'CAGE', 'CAME', 'CANE', 'CAPE', 'CARE', 'CASE',
    'CAVE', 'DAME', 'DANE', 'DARE', 'DATE', 'DAZE', 'FACE', 'FADE', 'FAKE',
    'FAME', 'FARE', 'FATE', 'GAVE', 'HALE', 'HAME', 'HANG', 'HARE', 'HATE',
    'HAVE', 'LAKE', 'LAME', 'LANE', 'LATE', 'MAZE', 'NAME', 'NAVE', 'PACE',
    'PACK', 'PALE', 'PANE', 'PANG', 'PARE', 'PAST', 'RAGE', 'RAKE', 'RANG',
    'RATE', 'RAVE', 'SAGE', 'SAKE', 'SALE', 'SAME', 'SANE', 'SAVE', 'TAKE',
    'TALE', 'TAME', 'TANG', 'TAPE', 'VALE', 'VANE', 'VANG', 'VARY', 'WAGE',
    'WAKE', 'WANE', 'WANG'
])

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in words and new_word != word:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    if start not in words or end not in words:
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

# Find the path from 'MAGE' to 'PANG'
path = find_path('MAGE', 'PANG')
if path:
    print(','.join(path))
else:
    print("No path found")