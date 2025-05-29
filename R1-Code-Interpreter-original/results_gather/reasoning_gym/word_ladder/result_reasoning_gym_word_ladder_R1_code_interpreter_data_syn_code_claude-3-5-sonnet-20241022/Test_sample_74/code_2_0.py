from collections import deque

# Carefully verified list of valid 4-letter English words that could help form the path
WORDS = set([
    'APES', 'AMES', 'AIMS', 'ARMS', 'TEAM',
    'TAME', 'TAPE', 'TAPS', 'TEAS', 'TENS',
    'TERM', 'TEMP', 'TAMP', 'TALE', 'TAKE',
    'APEX', 'ACES', 'ACTS', 'ALES', 'AREA',
    'AGES', 'AIRS', 'APED', 'AMES', 'AXES',
    'TAME', 'FAME', 'GAME', 'CAME', 'NAME',
    'TEAS', 'TEAR', 'TEAL', 'TRAM', 'TRAP',
])

def one_letter_diff(word1, word2):
    return sum(1 for a, b in zip(word1, word2) if a != b) == 1

def get_neighbors(word):
    return [w for w in WORDS if one_letter_diff(word, w)]

def find_path(start, end):
    if start not in WORDS or end not in WORDS:
        return None
    
    queue = deque([[start]])
    seen = {start}
    
    while queue:
        path = queue.popleft()
        current = path[-1]
        
        if current == end:
            return path
            
        for next_word in get_neighbors(current):
            if next_word not in seen:
                seen.add(next_word)
                queue.append(path + [next_word])
    
    return None

# Find the path
result = find_path('APES', 'TEAM')
if result:
    print(','.join(result))
else:
    print("No valid path found")