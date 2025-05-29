from collections import deque, defaultdict

# Define a set of common 4-letter words
WORDS = set([
    'APES', 'APEX', 'APES', 'APPS', 'APES', 'ALES', 'ACES', 'ACTS',
    'ARMS', 'ARMY', 'AREA', 'APEX', 'ACES', 'AGES', 'AIMS', 'AIRS',
    'TEAM', 'TEAR', 'TEAL', 'TRAM', 'TERM', 'TEMP', 'TAMP', 'TAPE',
    'TALE', 'TAKE', 'TAME', 'SAME', 'CAME', 'TAPS', 'TOPS', 'TOES',
    'TIES', 'TINS', 'TENS', 'TERN', 'TEEN', 'THAN', 'THAT', 'THEM',
    'APES', 'APEX', 'ACES', 'ACTS', 'ARMS', 'ALES', 'AREA', 'AGES'
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
    if start not in WORDS or end not in WORDS:
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

# Find the path from 'APES' to 'TEAM'
path = find_path('APES', 'TEAM')
if path:
    print(','.join(path))
else:
    print("No valid path found")