from collections import deque

# Carefully curated list focusing on common 4-letter words
word_list = '''SHED SHOE SHOD SHOP CHOP CHIP CHIT CHAT THAT THAN THEN WHEN WHET
WHAT CHAT COAT GOAT GNAT THAT THAN THIN CHIN CHIC CHIP CHOP SHOP SHOT SHUT SHUN
SHIN SKIN SKIP SHIP CHIP CHOP CROP CRAP CRIB SHIP SHOP CHOP CHIP CLIP CRIB SHED
SHELL SPELL SPILL STILL STALL STALE SHALE SHAKE SHARE SCARE SCORE SCORN SWORN
SHORN SHORE SHARE SHADE SPADE SPACE SPARE SCARE SCORE STORE STORK STARK STACK
SLACK SHACK SHARK SHARP HARP HARD CARD CORD CURD CURB CRIB'''

word_set = set(word_list.split())

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
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
        if word == end:
            return path
        
        for neighbor in get_neighbors(word):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    return None

# Find and print the path
path = find_path('SHED', 'CRIB')
if path:
    print(','.join(path))
else:
    print("No path found")