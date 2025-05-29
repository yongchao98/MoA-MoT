def get_word_list():
    # Larger set of common 4-letter words
    return set('''
    TIFF TOFF TUFF TAFT TUFT TILT TILL TILE TALE TAME FAME FACE FADE FATE GATE GAVE GAZE MAZE
    RIFF RIFT RAFT RATE MATE MAZE
    TIFF TIFT LIFT LIFE LIKE LAKE MAKE MAZE
    TIFF MIFF MIFF MILE MALE MAZE
    TIFF TIFT SIFT SOFT SORT SORE SARE SAFE SAME MAZE
    TIFF TIFT GIFT GIST GUST MUST MUSE MAZE
    TIFF TIFT LIFT LIST LOST LOSE LORE MARE MAZE
    TIFF TIFT RIFT RAFT RAVE RASE ROSE RAZE MAZE
    TIFF TIFT SIFT SITE SATE SAFE SAME MAZE
    '''.split())

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
from collections import deque
word_set = get_word_list()
path = find_path('TIFF', 'MAZE', word_set)

if path:
    print(','.join(path))
else:
    print("No path found")