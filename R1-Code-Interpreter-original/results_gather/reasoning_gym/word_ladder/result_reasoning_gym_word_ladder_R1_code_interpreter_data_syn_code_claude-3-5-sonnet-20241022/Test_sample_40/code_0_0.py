from collections import deque

# Common 4-letter words in English
words = set([
    'USED', 'FILE', 'USE', 'FUSE', 'FUEL', 'FILL', 'FIRE', 'FINE', 'FADE',
    'FACE', 'FATE', 'FAME', 'FILE', 'FIVE', 'FISH', 'FIST', 'FAST', 'FAIL',
    'FALL', 'FEEL', 'FEED', 'FLED', 'FLEE', 'FLEW', 'FLEX', 'FLEA', 'FOLD',
    'FOLK', 'FOOD', 'FOOL', 'FOOT', 'FORD', 'FORE', 'FORK', 'FORM', 'FORT',
    'FOUL', 'FOUR', 'FREE', 'FULL', 'FUND', 'FURY', 'FUSE', 'USEE', 'USER',
    'USES', 'USED', 'FILE'
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