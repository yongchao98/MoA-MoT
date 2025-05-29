from collections import deque, defaultdict

def get_word_list():
    # Common 4-letter words
    words = '''TIFF TILE TILL TALL TALE MALE MAZE
               TIME TINE MINE MAKE MACE
               TIFF RIFF RIFT RAFT RAVE MAZE
               LIFE LIKE MIKE MAKE MAZE'''.split()
    return set(words)

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            neighbors.append(new_word)
    return neighbors

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
                queue.append(path + [neighbor])
                seen.add(neighbor)
    
    return None

# Find the path
word_set = get_word_list()
path = find_path('TIFF', 'MAZE', word_set)

if path:
    print(','.join(path))
else:
    print("No path found")