from collections import deque, defaultdict

# Define a set of common 4-letter words
WORDS = set([
    'LOUD', 'LOAD', 'LEAD', 'LEND', 'LAND', 'LAID', 'LAID', 'TELL', 'TALL', 'TEAL',
    'TELL', 'TILL', 'TOLL', 'TOLD', 'LOAD', 'LORD', 'LORE', 'LOSE', 'LOST', 'LEST',
    'LAST', 'LATE', 'TALE', 'TALK', 'TALL', 'TELL', 'SELL', 'SEAL', 'TEAL', 'TEAM',
    'TEAR', 'TEND', 'TENT', 'TEST', 'TEXT', 'THAT', 'THAN', 'THEM', 'THEN', 'LEAD',
    'LEAN', 'LEAP', 'LEAR', 'LEND', 'LENS', 'LENT', 'LESS', 'LEST', 'LETS'
])

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in WORDS and new_word != word:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    if start == end:
        return [start]
    
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        word = path[-1]
        
        for neighbor in get_neighbors(word):
            if neighbor not in visited:
                visited.add(neighbor)
                new_path = list(path)
                new_path.append(neighbor)
                
                if neighbor == end:
                    return new_path
                    
                queue.append(new_path)
    
    return None

# Find the path and print it
path = find_path('LOUD', 'TELL')
if path:
    print(','.join(path))
else:
    print("No path found")