from collections import deque

# More comprehensive list of 4-letter English words
words = {
    'READ', 'ROAD', 'RAID', 'RAIL', 'SAIL', 'SHIP',
    'HEAD', 'HEAL', 'HEAR', 'HEAT', 'HEAP',
    'LEAD', 'LOAD', 'LEAN', 'LEAP',
    'REAL', 'SEAL', 'SEAT', 'SLIP', 'SHOP',
    'REED', 'SEED', 'SEND', 'SHED', 'SHAD',
    'REAP', 'HEAP', 'HELP', 'HELD',
    'READ', 'ROAD', 'ROAM', 'FOAM', 'FEAR',
    'BEAD', 'BEAM', 'SEAM', 'SLAM', 'SLAP',
    'READ', 'REAL', 'SEAL', 'SEAM', 'SLAM', 'SLAP', 'SHIP',  # Additional path
    'HEAD', 'HERD', 'HARD', 'HAND', 'SAND', 'SEND', 'SENT',
    'RIDE', 'SIDE', 'SITE', 'SIRE', 'HIRE', 'HIDE', 'HIKE',
    'REAR', 'BEAR', 'BEAT', 'BELT', 'BOLT', 'BOAT', 'COAT'
}

def get_neighbors(word):
    neighbors = set()
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in words and new_word != word:
                neighbors.add(new_word)
    return neighbors

def find_path(start, end):
    if start not in words or end not in words:
        return None
    
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        current = path[-1]
        
        if current == end:
            return path
        
        for next_word in get_neighbors(current):
            if next_word not in visited:
                visited.add(next_word)
                new_path = list(path)
                new_path.append(next_word)
                queue.append(new_path)
    
    return None

# Find and print the path
result = find_path('READ', 'SHIP')
if result:
    print(','.join(result))
else:
    print("No path found")