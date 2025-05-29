from collections import deque, defaultdict

# First, let's create a set of common 4-letter words
words = set('''
JADE JADE JAKE TAKE TALE TALL TELL TILL TILL TILL TOLL TOLE TOSS
BADE FADE GATE HATE LATE MATE PALE SALE TALE WADE
JACK PACK BACK BARK DARK DARE DATE DATE FAME GAME LAME NAME SAME TAME
TOLD FOLD HOLD SOLD COLD BOLD GOLD
TOSS LOSS LESS TESS TEAS TEAR TEAL TELL TALL TAIL TAIL
'''.split())

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in words and new_word != word:
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
path = find_path('JADE', 'TOSS')
if path:
    print(','.join(path))
else:
    print("No path found")