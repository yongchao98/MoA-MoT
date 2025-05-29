from collections import deque

# Expanded 4-letter words dictionary
words = set([
    'USED', 'FILE', 'FUSE', 'FUEL', 'FILL', 'FIRE', 'FINE', 'FADE',
    'FACE', 'FATE', 'FAME', 'FILE', 'FIVE', 'FISH', 'FIST', 'FAST', 'FAIL',
    'FALL', 'FEEL', 'FEED', 'FLED', 'FLEE', 'FLEW', 'FLEX', 'FLEA', 'FOLD',
    'FOLK', 'FOOD', 'FOOL', 'FOOT', 'FORD', 'FORE', 'FORK', 'FORM', 'FORT',
    'FOUL', 'FOUR', 'FREE', 'FULL', 'FUND', 'FURY', 'FUSE', 'USER', 'USES',
    'USED', 'FILE', 'FIDE', 'FIKE', 'FILA', 'FILO', 'FILS', 'FIND', 'FINK',
    'FISE', 'FUSS', 'FUSE', 'FUDE', 'FADE', 'FAME', 'FARE', 'FATE', 'FAVE',
    'FAZE', 'FEED', 'FEEL', 'FEES', 'FEET', 'FELL', 'FELT', 'FEND', 'FERN',
    'FETE', 'FEUD', 'FIED', 'FIFE', 'FILA', 'FILE', 'FILL', 'FILM', 'FIND',
    'FINE', 'FIRE', 'FIRM', 'FIST', 'FIVE', 'FIXE', 'FIZE', 'FLEA', 'FLED',
    'FLEE', 'FLEW', 'FLEX', 'FLEY', 'FLIE', 'FLIP', 'FLIR', 'FOLD', 'FOLK',
    'FOND', 'FONT', 'FOOD', 'FOOL', 'FOOT', 'FORD', 'FORE', 'FORK', 'FORM',
    'FORT', 'FOUL', 'FOUR', 'FREE', 'FRET', 'FROE', 'FROM', 'FUEL', 'FULL',
    'FUME', 'FUND', 'FUNK', 'FURY', 'FUSE', 'FUSS', 'USEE', 'USER', 'USES',
    'USE', 'UFED', 'URED', 'UDED', 'UDES', 'UDOS', 'URED', 'URGE', 'URNS',
    'USED'
])

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in words and new_word != word:
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

# Find the path
path = find_path('USED', 'FILE')
if path:
    print(','.join(path))
else:
    print("No solution found")