from collections import deque

# Expanded list of common 4-letter words that might help form the path
words = '''
GHAT THAT CHAT WHAT GOAT COAT BOAT BEAT BELT BOLT BULK BUCK TUCK
THAT THAN THEN THIN THIS THUS TUCK TICK TACK TALK TANK TUNK
GHAT CHAT THAT THAN THEN THEM TICK TUCK TACK
GHAT THAT TEAT TENT TEXT TUCK TACK TALK TANK
BACK BUCK BANK TANK TALK WALK HAWK HARK HARD HAND BAND BANK
THAT WHAT CHAT COAT GOAT GHAT
TICK TICK TACK TALK TANK TUNK TUCK
'''.split()

def get_valid_words(length=4):
    return set(word for word in words if len(word) == length)

def get_neighbors(word, valid_words):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in valid_words:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    valid_words = get_valid_words()
    if start not in valid_words or end not in valid_words:
        return None
    
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        word = path[-1]
        
        if word == end:
            return path
            
        for neighbor in get_neighbors(word, valid_words):
            if neighbor not in visited:
                visited.add(neighbor)
                new_path = list(path)
                new_path.append(neighbor)
                queue.append(new_path)
    
    return None

# Find and print the path
path = find_path('GHAT', 'TUCK')
if path:
    print(','.join(path))
else:
    print("No valid path found")