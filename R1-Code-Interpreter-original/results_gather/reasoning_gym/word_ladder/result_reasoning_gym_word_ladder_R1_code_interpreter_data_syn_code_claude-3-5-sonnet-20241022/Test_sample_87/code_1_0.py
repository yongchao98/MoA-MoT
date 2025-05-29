from collections import deque

# Define a more comprehensive set of 4-letter words that we know exist in English
words = {
    'FAVS', 'FADS', 'SUDS', 'SUBS', 'SUES', 'SAVE', 'FACE', 'FADE', 'FAME', 'FARE',
    'FAST', 'FATE', 'FATS', 'HAVE', 'CAVE', 'DAVE', 'FAKE', 'FALL', 'FAME', 'FARE',
    'FARM', 'FAST', 'FATE', 'FAVE', 'FAWN', 'FEAR', 'FEAT', 'FEED', 'FEEL', 'FEES',
    'FEET', 'FELL', 'FELT', 'FEND', 'FERN', 'FEST', 'FEUD', 'FADS', 'SUDS', 'SADS',
    'SAYS', 'SEAS', 'SEAM', 'SEAL', 'SEEP', 'SEES', 'SEND', 'SENT', 'SETS', 'SEWS',
    'SUBS', 'SUCH', 'SUCK', 'SUED', 'SUES', 'SUIT', 'SUMS', 'SUNG', 'SUNK', 'SUNS',
    'SURE', 'SURF', 'SUSS', 'FAUS', 'FAYS', 'SADS'
}

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in words:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    if start not in words or end not in words:
        return None
    
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
path = find_path('FAVS', 'SUDS')
if path:
    print(','.join(path))
else:
    print("No path found")