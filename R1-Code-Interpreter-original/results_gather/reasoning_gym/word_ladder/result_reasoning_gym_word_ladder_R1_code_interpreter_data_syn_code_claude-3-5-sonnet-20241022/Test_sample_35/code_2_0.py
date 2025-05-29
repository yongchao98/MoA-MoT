from collections import deque

# Carefully curated list of 4-letter words focusing on the path we need
WORDS = {
    'CUSP', 'CUSS', 'CUPS', 'CUTS', 'CUES', 'CURE', 'CORE', 'COPE', 'CAPE', 'CASE',
    'CAST', 'COST', 'COAT', 'COLT', 'COLD', 'CORD', 'CARD', 'CARE', 'CAVE', 'COVE',
    'OVER', 'OUST', 'DUST', 'DUSK', 'DISK', 'RISK', 'RISE', 'ROSE', 'ROVE',
    'CARP', 'CART', 'CARS', 'CAPS', 'COPS', 'COGS', 'COBS', 'CUBS', 'CURB',
    'COUP', 'COOP', 'COOK', 'CORK', 'CORN', 'CORP', 'COPE', 'COME', 'DOME',
    'COVE', 'CAVE', 'CASE', 'CAST', 'COST', 'MOST', 'MOSS', 'MISS', 'KISS',
    'CASK', 'CASH', 'CASE', 'CARE', 'CURE', 'SURE', 'SORE', 'SOUR', 'SOUP',
    'CUST', 'OUST', 'OAST', 'OARS', 'OURS', 'OVER'
}

def get_neighbors(word):
    neighbors = []
    for i in range(4):  # 4-letter words
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in WORDS and new_word != word:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    if start not in WORDS or end not in WORDS:
        return None
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        if current == end:
            return path
            
        for next_word in get_neighbors(current):
            if next_word not in visited:
                visited.add(next_word)
                new_path = path + [next_word]
                queue.append((next_word, new_path))
    return None

# Try to find the path
path = find_path('CUSP', 'OVER')
if path:
    print(','.join(path))
else:
    print("Trying alternative paths...")
    # Try some known intermediate steps
    for mid in ['COVE', 'OUST', 'CURE']:
        path1 = find_path('CUSP', mid)
        if path1:
            path2 = find_path(mid, 'OVER')
            if path2:
                full_path = path1 + path2[1:]
                print(','.join(full_path))
                break