from collections import deque, defaultdict

# Common 4-letter words
words = '''SCAB SLAB SLOB SLOT SLIT SLIP SLID SIDE
          STAB STAR SCAR SCAN SPAN SPIN SPIT SLED
          SEED SAID SAND SEND SEAT BEAT BELT BENT
          SHED SHOE SHOP SHOT SPOT SOOT SOFT SORT'''.split()

# Create adjacency list
def get_neighbors(word):
    neighbors = []
    for w in words:
        diff = sum(1 for a, b in zip(word, w) if a != b)
        if diff == 1:
            neighbors.append(w)
    return neighbors

# BFS to find shortest path
def find_path(start, end):
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        word = path[-1]
        
        if word == end:
            return ','.join(path)
            
        for neighbor in get_neighbors(word):
            if neighbor not in visited:
                visited.add(neighbor)
                new_path = list(path)
                new_path.append(neighbor)
                queue.append(new_path)
    
    return None

result = find_path('SCAB', 'SIDE')
print(result)