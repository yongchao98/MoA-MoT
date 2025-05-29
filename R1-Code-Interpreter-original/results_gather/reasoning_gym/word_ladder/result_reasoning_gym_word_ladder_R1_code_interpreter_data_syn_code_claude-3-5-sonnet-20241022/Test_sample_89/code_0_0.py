from collections import deque, defaultdict

# Common 4-letter words in English
words = '''GHAT THAT CHAT WHAT GOAT COAT BOAT BEAT BELT BOLT BULK BUCK TUCK
          THAT THAN THEN THIN THIS THUS TUCK
          GHAT CHAT THAT THAN THEN THEM TICK TUCK
          GHAT THAT TEAT TENT TENT TENT TEXT TUCK'''.split()

# Create set for O(1) lookup
word_set = set(words)

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in word_set:
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
path = find_path('GHAT', 'TUCK')
if path:
    print(','.join(path))
else:
    print("No valid path found")