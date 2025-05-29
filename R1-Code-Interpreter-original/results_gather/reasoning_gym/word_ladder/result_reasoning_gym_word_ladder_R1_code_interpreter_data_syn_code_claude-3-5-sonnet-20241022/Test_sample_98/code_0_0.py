from collections import deque
import urllib.request

# Download word list
word_url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
response = urllib.request.urlopen(word_url)
word_list = set(word.decode('utf-8').strip().upper() for word in response.readlines())

# Filter only 4-letter words
valid_words = set(word for word in word_list if len(word) == 4)

def get_neighbors(word):
    neighbors = []
    for i in range(4):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in valid_words:
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

# Find the path
path = find_path('HYPE', 'PUCK')
if path:
    print(','.join(path))
else:
    print("No valid path found")