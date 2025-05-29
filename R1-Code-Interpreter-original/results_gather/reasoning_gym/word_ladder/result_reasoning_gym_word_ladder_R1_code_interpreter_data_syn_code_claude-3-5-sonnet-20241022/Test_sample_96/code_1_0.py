from collections import deque

# Expanded set of 4-letter words
WORDS = set([
    'OPUS', 'ONUS', 'OUTS', 'OPTS', 'OATS', 'OILS', 'OWLS', 'OWNS', 'OAKS',
    'SKIP', 'SKIN', 'SKID', 'SKIS', 'SKIT', 'SLID', 'SLIM', 'SLIT', 'SLIP',
    'SPIN', 'SPIT', 'SPOT', 'SPAT', 'SPAN', 'SPAR', 'SPAS', 'SPED',
    'THUS', 'THIS', 'THIN', 'THAN', 'THAT', 'THEM', 'THEY',
    'PLUS', 'PUSS', 'PUPS', 'PUTS', 'PITS', 'PINS', 'PENS', 'PETS',
    'SOUS', 'SOUR', 'SOUP', 'SOUL', 'SOLS', 'SOIL', 'SOAR', 'SOAK',
    'SODS', 'SOWS', 'SONS', 'SONG', 'SORT', 'SORE', 'SOLE', 'SOLD',
    'SITS', 'SITE', 'SIRE', 'SINE', 'SIDE', 'SIDS', 'SICS', 'SICK',
    'SIPS', 'SIRS', 'SINS', 'SINK', 'SING', 'SINE', 'SIGN', 'SIGH',
    'OPUS', 'ONUS', 'OUTS', 'OPTS', 'OATS', 'OILS', 'OWLS', 'OWNS',
    'SKIS', 'SKIT', 'SKIP', 'SKIN', 'SKIM', 'SKID', 'SKEW', 'SKED',
    'SOIS', 'SODS', 'SOLS', 'SOUS', 'SOWS', 'SONS', 'SOTS', 'SOHS'
])

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
                new_path = list(path)
                new_path.append(neighbor)
                queue.append(new_path)
    
    return None

# Find the path and print the result
path = find_path('OPUS', 'SKIS')
if path:
    print(','.join(path))
else:
    print("No valid path found")