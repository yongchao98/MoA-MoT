from collections import deque, defaultdict

# Create a set of valid 4-letter English words
words = set()
try:
    with open('words_alpha.txt', 'r') as f:
        for line in f:
            word = line.strip().upper()
            if len(word) == 4:
                words.add(word)
except:
    # If file not available, use a minimal set of common 4-letter words
    words = {'MODE', 'MOVE', 'MOLE', 'MALE', 'MARE', 'MORE', 'MOTE', 'MADE', 'CODE', 'CORE',
             'TORE', 'TIRE', 'TYRE', 'TOME', 'TONE', 'TUNE', 'TALE', 'TAME', 'TAKE', 'TARE'}

def get_neighbors(word):
    neighbors = []
    for i in range(4):
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

# Find the path
path = find_path('MODE', 'TYRE')
if path:
    print(','.join(path))
else:
    print("No valid path found")