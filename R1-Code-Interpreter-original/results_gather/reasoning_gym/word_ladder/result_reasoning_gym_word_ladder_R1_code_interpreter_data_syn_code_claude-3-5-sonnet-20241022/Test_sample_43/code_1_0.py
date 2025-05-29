def can_transform(word1, word2):
    # Check if two words differ by exactly one letter
    diff = 0
    for c1, c2 in zip(word1, word2):
        if c1 != c2:
            diff += 1
    return diff == 1

# Define a more targeted set of common 4-letter words
words = '''TUMS TUBS TABS TAPS TIPS TOPS TOPE TAPE TALE PALE POLE POPE PIPE PREP
PENS PINS TINS TONS TOMS TOME TIME TIRE TARE PARE PURE PORE PYRE PREY
TAME DAME DARE PARE PORE PERK PEAK PEAL TEAL TEAM TEMP TAMP TRAP
PUMP PUMP PULP PUP PUPS PEPS PEEP DEEP DEER PEER PEAR TEAR TERM PERM'''.split()

def find_path(start, end, word_list):
    from collections import deque
    
    # Create set for O(1) lookup
    word_set = set(word_list)
    if start not in word_set or end not in word_set:
        return None
        
    # BFS queue and visited set
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        current = path[-1]
        
        # Check all possible neighbors
        for word in word_set:
            if word not in visited and can_transform(current, word):
                if word == end:
                    return path + [word]
                visited.add(word)
                queue.append(path + [word])
    
    return None

# Find path
path = find_path('TUMS', 'PREP', words)
if path:
    print(','.join(path))
else:
    print("No path found")