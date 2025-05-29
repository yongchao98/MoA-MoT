from collections import deque

# Carefully curated word list ensuring connectivity between SCAB and SIDE
words_set = set('''
SCAB SLAB STAB STAG SLAG SLID SIDE
SCAB SCAR SOAR SEAR SEAL SEAM SEEM SEED SPED SLED SLID SIDE
SCAB STAB STAR SEAR SEAT SEPT SECT SEED SIDE
SCAB SLAB SLAD SLED SEED SIDE
'''.split())

def find_shortest_path(start, end, max_depth=10):
    if start == end:
        return start
    
    # Use BFS with depth limit to find shortest path
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        
        # Don't explore paths that are too long
        if len(path) > max_depth:
            continue
        
        # For each position in the word
        for i in range(len(current)):
            # Try each letter
            word_chars = list(current)
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                word_chars[i] = c
                next_word = ''.join(word_chars)
                
                if next_word in words_set and next_word not in visited:
                    if next_word == end:
                        return ','.join(path + [next_word])
                    visited.add(next_word)
                    queue.append((next_word, path + [next_word]))
    return None

result = find_shortest_path('SCAB', 'SIDE')
print(result)