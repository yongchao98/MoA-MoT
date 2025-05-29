from collections import deque, defaultdict

# Create a set of 4-letter words
def get_word_set():
    words = set()
    try:
        with open('words_alpha.txt', 'r') as f:
            for line in f:
                word = line.strip().upper()
                if len(word) == 4:
                    words.add(word)
    except:
        # Fallback to a basic set of common 4-letter words
        words = {'FAVS', 'FADS', 'SUDS', 'SAVE', 'FACE', 'FADE', 'FAME', 'FARE',
                'FAST', 'FATE', 'FATS', 'HAVE', 'CAVE', 'DAVE', 'FADE', 'FAKE',
                'FALL', 'FAME', 'FARE', 'FARM', 'FAST', 'FATE', 'FAVE', 'FAWN',
                'FEAR', 'FEAT', 'FEED', 'FEEL', 'FEES', 'FEET', 'FELL', 'FELT',
                'FEND', 'FERN', 'FEST', 'FEUD', 'FADS', 'SUDS'}
    return words

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            neighbors.append(new_word)
    return neighbors

def find_path(start, end, words):
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
            if neighbor in words and neighbor not in visited:
                visited.add(neighbor)
                new_path = list(path)
                new_path.append(neighbor)
                queue.append(new_path)
    
    return None

# Find the path
words = get_word_set()
path = find_path('FAVS', 'SUDS', words)

if path:
    print(','.join(path))
else:
    print("No path found")