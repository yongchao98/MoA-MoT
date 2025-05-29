from collections import deque

# Using a more carefully selected set of common 4-letter words
words = set('''
REIN RAIN GAIN GRIN GRIT GRID GRIS GRIM GRIM GRIM TRIM THIN CHIN CHIS
REIN RUIN RUNT RANT RANT CANT CHAT CHIT CHIN CHIS
REIN RAIN RAIS RAYS RATS CATS CHATS CHITS CHIS
REIN REIN VEIN VAIN VAIN CHIN CHIS
'''.split())

def hamming_distance(word1, word2):
    return sum(1 for a, b in zip(word1, word2) if a != b)

def get_neighbors(word):
    return [w for w in words if hamming_distance(word, w) == 1]

def find_path(start, end):
    if start not in words or end not in words:
        return None
        
    queue = deque([(start, [start])])
    seen = {start}
    
    while queue:
        current, path = queue.popleft()
        if current == end:
            return path
            
        for next_word in get_neighbors(current):
            if next_word not in seen:
                seen.add(next_word)
                queue.append((next_word, path + [next_word]))
    return None

# Find the shortest path
result = find_path('REIN', 'CHIS')
if result:
    print(','.join(result))
else:
    print("No valid path found")