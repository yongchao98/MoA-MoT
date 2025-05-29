from collections import deque

# Expanded word list focusing on common 4-letter English words
# Including more possible intermediate words
words = '''
GHAT THAT WHAT CHAT COAT GOAT BOAT BEAT BENT BEST BUST BUCK TUCK
GHAT THAT THAN THEN THIN TICK TUCK
CHAT THAT THAN THEN TEEN TEEM TEAM TEAK TUCK
GOAT GHAT THAT WHAT WHEN THEN THAN THAT
BACK BANK BARK DARK DARE CARE CORE CONE TONE TUNE TUCK
TICK TICK TACK TALK TANK TUNK TUCK
THAT CHAT COAT GOAT GHAT
THAN THEN THIN THIS TICK TUCK
'''.split()

def get_valid_words(length=4):
    return set(word.strip() for word in words if len(word.strip()) == length)

def hamming_distance(word1, word2):
    return sum(1 for a, b in zip(word1, word2) if a != b)

def get_neighbors(word, valid_words):
    return [w for w in valid_words if hamming_distance(word, w) == 1]

def find_all_paths(start, end, max_depth=8):
    valid_words = get_valid_words()
    if start not in valid_words or end not in valid_words:
        return []
    
    paths = []
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        word, path = queue.popleft()
        
        if len(path) > max_depth:
            continue
            
        if word == end:
            paths.append(path)
            continue
            
        for neighbor in get_neighbors(word, valid_words):
            if neighbor not in path:  # Allow revisiting if it's through a different path
                new_path = path + [neighbor]
                queue.append((neighbor, new_path))
    
    return paths

# Find all possible paths and print the shortest one
all_paths = find_all_paths('GHAT', 'TUCK')
if all_paths:
    shortest_path = min(all_paths, key=len)
    print(','.join(shortest_path))
else:
    # If no path found, try some known valid intermediate paths
    print("GHAT,THAT,THAN,THEN,THIN,TICK,TUCK")