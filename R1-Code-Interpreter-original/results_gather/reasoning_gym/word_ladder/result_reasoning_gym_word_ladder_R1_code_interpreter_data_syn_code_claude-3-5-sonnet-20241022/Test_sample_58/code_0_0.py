from collections import deque, defaultdict

# Common 4-letter English words
words = '''CLAN CLAD CLAY SLAY SLAP SNAP SOAP SOUP SOUR RUSH
          PLAN PLAY SLAT FLAT FLAP CLAP CRAP TRAP TRIP SLIP
          SLAM CLAW FLAW FLAN SCAN SPAN SPIN SPIT SPUR SPUN
          CRAB GRAB GRAM GRIM GRIT GRID RAID RAIN RUIN RASH
          CASH BASH BUSH RUSE FUSE RUSK'''.split()

def get_neighbors(word):
    neighbors = []
    for i in range(4):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in word_set and new_word != word:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    if start == end:
        return [start]
    
    word_set.add(start)
    word_set.add(end)
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        word, path = queue.popleft()
        for neighbor in get_neighbors(word):
            if neighbor == end:
                return path + [end]
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    return None

# Convert words list to set for O(1) lookup
word_set = set(words)

# Find the path
path = find_path('CLAN', 'RUSH')

if path:
    print(','.join(path))
else:
    print("No path found")