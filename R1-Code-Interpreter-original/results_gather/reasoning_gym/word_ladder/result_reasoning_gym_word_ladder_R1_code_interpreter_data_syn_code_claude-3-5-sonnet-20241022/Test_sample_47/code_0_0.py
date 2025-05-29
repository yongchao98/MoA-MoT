from collections import deque, defaultdict

# Create a set of valid 4-letter words
valid_words = set()
try:
    with open('/usr/share/dict/words', 'r') as f:
        for line in f:
            word = line.strip().upper()
            if len(word) == 4 and word.isalpha():
                valid_words.add(word)
except:
    # Fallback to a basic set of common 4-letter words
    valid_words = {'PROD', 'PROM', 'PRAM', 'MRAM', 'MARL',
                  'PARD', 'PARK', 'MARK', 'PROM', 'FROM',
                  'PRAM', 'GRAM', 'PROM', 'PROP', 'PROW',
                  'PRAY', 'PRAT', 'PLOD', 'PLOT', 'PROT'}

def get_neighbors(word):
    neighbors = []
    for i in range(4):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in valid_words:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    if start not in valid_words or end not in valid_words:
        return None
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        word, path = queue.popleft()
        if word == end:
            return path
            
        for neighbor in get_neighbors(word):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    
    return None

# Find the path and print it
path = find_path('PROD', 'MARL')
if path:
    print(','.join(path))
else:
    print("No valid path found")