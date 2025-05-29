from collections import deque

# Carefully selected set of 4-letter English words focusing on possible paths
WORDS = {
    'FEED', 'REED', 'NEED', 'SEED',
    'FEEL', 'PEEL', 'REEL', 'HEEL',
    'FELL', 'TELL', 'SELL', 'WELL',
    'FLAN', 'PLAN', 'CLAN', 'ELAN',
    'LEAN', 'LEAD', 'LEAK', 'LEAP',
    'FLED', 'FEED', 'FEES', 'FEET',
    'FOLD', 'FOND', 'FIND', 'FINE',
    'FELD', 'HELD', 'HELP', 'HEAP',
    'HEAD', 'HEAL', 'HEAR', 'HEAT',
    'FEAT', 'FEAR', 'DEAR', 'DEAN'
}

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in WORDS and new_word != word:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        word = path[-1]
        
        if word == end:
            return path
            
        for neighbor in get_neighbors(word):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(path + [neighbor])
    
    return None

# Find and print the path
result = find_path('FEED', 'ELAN')
if result:
    print(','.join(result))
else:
    print("No valid path found")