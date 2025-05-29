from collections import deque, defaultdict

# Define a function to get valid 4-letter words
def get_word_list():
    import urllib.request
    url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
    response = urllib.request.urlopen(url)
    words = response.read().decode().split()
    return set(word.upper() for word in words if len(word) == 4)

# Function to get words that differ by one letter
def get_neighbors(word):
    neighbors = []
    for i in range(4):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            neighbors.append(new_word)
    return neighbors

# BFS to find the path
def find_path(start, end, word_set):
    queue = deque([[start]])
    seen = {start}
    
    while queue:
        path = queue.popleft()
        word = path[-1]
        
        if word == end:
            return path
            
        for neighbor in get_neighbors(word):
            if neighbor in word_set and neighbor not in seen:
                seen.add(neighbor)
                new_path = list(path)
                new_path.append(neighbor)
                queue.append(new_path)
    
    return None

# Main execution
words = get_word_list()
path = find_path('DAMS', 'DIED', words)

if path:
    print(','.join(path))
else:
    print("No solution found")