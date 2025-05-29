from collections import deque

# A carefully curated list of common 4-letter English words
# focusing on possible paths between APES and TEAM
WORDS = {
    'APES', 'APEX', 'ACES', 'AGES', 'ALES', 'AMES', 'ARMS',
    'TAPS', 'TAPE', 'TAME', 'TEAM', 'TEMP', 'TERM',
    'TEAS', 'TENS', 'TERN', 'THEM', 'TIES', 'TIPS',
    'TOPS', 'TOES', 'TOMS', 'TONS', 'TANS', 'TABS',
    'ARES', 'ARTS', 'AIMS', 'AXIS', 'AXES', 'ANTS'
}

def get_neighbors(word):
    neighbors = []
    for i in range(4):  # for each position
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in WORDS and new_word != word:
                neighbors.append(new_word)
    return neighbors

def find_shortest_path(start, target):
    if start not in WORDS or target not in WORDS:
        return None
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current_word, path = queue.popleft()
        
        if current_word == target:
            return path
        
        for neighbor in get_neighbors(current_word):
            if neighbor not in visited:
                visited.add(neighbor)
                new_path = path + [neighbor]
                queue.append((neighbor, new_path))
    
    return None

# Find the path
path = find_shortest_path('APES', 'TEAM')
if path:
    print(','.join(path))
else:
    print("No valid path found")