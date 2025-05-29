from collections import deque
import urllib.request

# Download word list
url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
response = urllib.request.urlopen(url)
words = set(word.decode('utf-8').strip().upper() for word in response if len(word.decode('utf-8').strip()) == 4)

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in words:
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
path = find_path('LAZE', 'NUKE')
if path:
    print(','.join(path))
else:
    print("No path found")