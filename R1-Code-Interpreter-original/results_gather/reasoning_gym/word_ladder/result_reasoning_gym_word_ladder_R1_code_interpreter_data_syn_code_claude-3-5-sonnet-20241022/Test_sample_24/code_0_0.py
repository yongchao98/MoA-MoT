from collections import deque, defaultdict
import urllib.request

# Download word list
url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
response = urllib.request.urlopen(url)
words = set(word.decode('utf-8').strip().upper() for word in response if len(word.strip()) == 4)

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in words and new_word != word:
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
path = find_path('ROPY', 'SODS')
if path:
    print(','.join(path))
else:
    print("No path found")