from collections import deque, defaultdict

# Define a set of common 4-letter English words
WORDS = {
    'BEAD', 'BEAN', 'BEAT', 'BEND', 'BEAR', 'BRED', 'DEAD', 'HEAD', 'LEAD', 'MEAD', 'READ',
    'BEAM', 'SEAM', 'SLAM', 'SPAN', 'SWAN', 'SAWN', 'DAWN', 'LAWN', 'PAWN', 'SAND', 'SEND',
    'SEAT', 'SEAN', 'SEEN', 'SEWN', 'SHED', 'SPED', 'SPIN', 'SAWN', 'SAID', 'SEAL', 'SEAM',
    'SEAR', 'SEED', 'SEEK', 'SEEM', 'SEEN', 'SEER', 'SEAS', 'SEAT', 'SEED', 'SEEK', 'SEEL'
}

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in WORDS:
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

# Find the path and print it
path = find_path('BEAD', 'SAWN')
if path:
    print(','.join(path))
else:
    print("No path found")