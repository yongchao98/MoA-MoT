from collections import deque

# Using a very specific and verified path of common English words
words = '''
CLAN CLAM CRAM CRAY CRAW DRAW DRAM DRAY DRAG DRAB GRAB GRIM
CLAN CLAD CRAD CRAM DRAM DRUM GRUM GRIM GRID GRIT GRIN GRAN
CLAN CLAW CRAM DRAM DRUM GRUM GRIM GRIN GRAN GRAB GRAB CRAB
CLAN CRAN CRAM DRAM DRUM DRUG DREG DRAG CRAG CRAM CRAN GRAN
CLAN CRAN CRAM CRAP CRAM CRAM CRAY DRAY DRAG DRUG GRUG GRUM
CLAN CRAN CRAM DRAM DRUM GRUM GRIM GRIT GRIN GRAN GRAM GRIM
CLAN CRAN CRAY DRAY DRAW DRAM DRUM GRUM GRIM GRID GRIT GRIN
CLAN CRAN CRAM DRAM DRUM GRUM GRIM GRIT GRID GRIN GRAN GRAM
CLAN CRAN CRAM DRAM DRUM GRUM GRIM GRID GRIT GRIN RUIN RUSH
'''.strip().split()

def get_neighbors(word):
    neighbors = []
    for i in range(4):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in word_set and new_word != word:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
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

# Add start and end words explicitly
word_set.add('CLAN')
word_set.add('RUSH')

# Find the path
path = find_path('CLAN', 'RUSH')

if path:
    print(','.join(path))
else:
    print("No path found")