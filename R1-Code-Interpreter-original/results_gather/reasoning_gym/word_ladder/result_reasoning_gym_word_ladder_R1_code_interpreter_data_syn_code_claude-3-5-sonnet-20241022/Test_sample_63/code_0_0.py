from collections import deque, defaultdict

# Define a set of common 4-letter English words
WORDS = {
    'SEED', 'SEND', 'SAND', 'BAND', 'BADE', 'BEAD', 'BEND', 'BEAT',
    'BATE', 'BARE', 'BASE', 'CASE', 'FADE', 'MADE', 'WADE', 'LADE',
    'LEAD', 'READ', 'ROAD', 'LOAD', 'DEED', 'FEED', 'NEED', 'WEED',
    'BEEN', 'SEEN', 'BEAN', 'MEAN', 'LEAN', 'DEAD', 'HEAD', 'BALD',
    'BOLD', 'SOLD', 'SAID', 'PAID', 'LAID', 'LAND', 'HAND', 'HARD',
    'CARD', 'CARE', 'DARE', 'DATE', 'GATE', 'HATE', 'LATE', 'MATE',
    'RATE', 'SAVE', 'SAME', 'SOME', 'COME', 'NAME', 'GAME', 'GAVE'
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
    if start not in WORDS or end not in WORDS:
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

# Find the path and print the result
path = find_path('SEED', 'BADE')
if path:
    print(','.join(path))
else:
    print("No valid path found")