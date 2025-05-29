from collections import deque, defaultdict

# Define a set of common 4-letter words
words = set([
    'WOLF', 'WOLD', 'WORD', 'WARD', 'WARN', 'WARN', 'GRAN', 'GRAM', 'GOLD', 'GOLF',
    'WOOL', 'WOOD', 'GOOD', 'GOAD', 'GRAD', 'GRAM', 'GRIM', 'GRID', 'GAIN', 'RAIN',
    'WOOF', 'WORE', 'GORE', 'GONE', 'GENE', 'GENT', 'WENT', 'WANT', 'WAND', 'WIND',
    'WING', 'WINE', 'WIPE', 'WIRE', 'WIRY', 'WORN', 'WORK', 'WALK', 'WALL', 'WELL'
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

# Find the path from WOLF to GRAN
path = find_path('WOLF', 'GRAN')
if path:
    print(','.join(path))
else:
    print("No path found")