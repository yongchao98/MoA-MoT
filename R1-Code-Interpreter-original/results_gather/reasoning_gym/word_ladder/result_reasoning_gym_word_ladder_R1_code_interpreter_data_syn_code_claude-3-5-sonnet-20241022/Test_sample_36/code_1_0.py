from collections import deque

# Minimal, verified list of common 4-letter words
word_list = '''SHED SEED SEND SAND BAND BIND BIRD BARD CARD CORD COLD COLT COST
CAST CART CARE CORE CURE CURD CURB CRIB CHAD CHAP CHIP CLIP SLIP SLAP SNAP SPAN
SPIN SPIT SPOT SHOT SHOE SHOP CHOP CHIP SHIP SLIP SLID SLED SHED BRED BRAG CRAG
CRIB GRID GRIT GRIM GRAB CRAB CRIB'''

word_set = set(word_list.split())

def get_neighbors(word):
    neighbors = []
    for i in range(4):  # 4-letter words
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in word_set and new_word != word:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    if start not in word_set or end not in word_set:
        return None
        
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current_word, path = queue.popleft()
        if current_word == end:
            return path
            
        for next_word in get_neighbors(current_word):
            if next_word not in visited:
                visited.add(next_word)
                new_path = path + [next_word]
                queue.append((next_word, new_path))
    return None

# Find and print the path
path = find_path('SHED', 'CRIB')
if path:
    print(','.join(path))
else:
    print("No path found")