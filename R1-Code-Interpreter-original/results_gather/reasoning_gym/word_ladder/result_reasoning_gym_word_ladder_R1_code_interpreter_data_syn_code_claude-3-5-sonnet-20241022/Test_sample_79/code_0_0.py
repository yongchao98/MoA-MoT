from collections import deque, defaultdict

# Define a set of common 4-letter English words
WORDS = {
    'LEAP', 'LEAD', 'LOAD', 'LEND', 'LAND', 'LEAN', 'LEAS', 'SEED',
    'SEND', 'SAID', 'SEAL', 'SEAT', 'SEAS', 'SEAR', 'SOLD', 'SOLE',
    'SOME', 'SAME', 'SAVE', 'SEEN', 'SEEM', 'SEEP', 'SEES', 'SOAP',
    'HEAD', 'HEAP', 'HEAR', 'HEAT', 'HELP', 'HELD', 'HEAL', 'READ',
    'REAL', 'REAP', 'REAR', 'ROAD', 'ROAR', 'SLED', 'SLID', 'SLIP',
    'SLIT', 'SLOT', 'SLOW', 'SLAP', 'SLOP', 'SHOP', 'SHOW', 'SHOT',
    'SHOE', 'SHED', 'SEED'
}

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

# Find and print the path
path = find_path('LEAP', 'SEED')
if path:
    print(','.join(path))
else:
    print("No path found")