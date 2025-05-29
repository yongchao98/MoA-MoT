from collections import deque

# More focused set of common 4-letter words
words_set = set('''
SCAB SCAR SCAN SEAT SAID SIDE
SCAT SEAT SAID SIDE
SCAD SEAD SEED SAID SIDE
'''.split())

def find_shortest_path(start, end):
    queue = deque([(start, [start])])
    visited = {start}
    shortest_path = None
    shortest_length = float('inf')
    
    while queue:
        current, path = queue.popleft()  # Use popleft() for BFS to ensure shortest path
        
        # If we already found a shorter path, skip longer paths
        if len(path) >= shortest_length:
            continue
            
        if current == end:
            if len(path) < shortest_length:
                shortest_path = path
                shortest_length = len(path)
            continue
            
        # Try changing each position
        for i in range(len(current)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                next_word = current[:i] + c + current[i+1:]
                if next_word in words_set and next_word not in visited:
                    visited.add(next_word)
                    queue.append((next_word, path + [next_word]))
    
    return ','.join(shortest_path) if shortest_path else None

result = find_shortest_path('SCAB', 'SIDE')
print(result)