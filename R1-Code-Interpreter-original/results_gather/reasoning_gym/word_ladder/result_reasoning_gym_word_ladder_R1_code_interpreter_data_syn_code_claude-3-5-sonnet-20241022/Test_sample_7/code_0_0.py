from collections import deque, defaultdict
import urllib.request

# Download word list
word_list_url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
response = urllib.request.urlopen(word_list_url)
words = set(word.decode('utf-8').strip().upper() for word in response.readlines())

# Filter only 4-letter words
words_4 = set(word for word in words if len(word) == 4)

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in words_4:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    if start not in words_4 or end not in words_4:
        return None
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        word, path = queue.popleft()
        if word == end:
            return path
            
        for neighbor in get_neighbors(word):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    
    return None

# Find the path
path = find_path('HIYA', 'NEED')
if path:
    print(','.join(path))
else:
    print("No valid path found")