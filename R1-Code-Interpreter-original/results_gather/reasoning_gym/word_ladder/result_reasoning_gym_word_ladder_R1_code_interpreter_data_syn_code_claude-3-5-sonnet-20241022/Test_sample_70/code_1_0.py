from collections import deque

# Expanded set of 4-letter words
WORDS = set([
    'WEAN', 'MEAN', 'MEAT', 'MESA', 'MEAL', 'MEAD', 'MEAS', 'WEAR', 'WEAK',
    'WEED', 'WEEK', 'WEEP', 'WEST', 'MELT', 'MEND', 'MEET', 'MEEK', 'MEED',
    'MOAN', 'MEAS', 'META', 'METE', 'MEWL', 'SEAL', 'SEAM', 'SEAN', 'SEAT',
    'WELD', 'WELL', 'WELT', 'WENT', 'WERE', 'WEST', 'BEAN', 'LEAN', 'DEAN',
    'MEAD', 'MELD', 'MEND', 'MENS', 'MENT', 'MESS', 'MOST', 'MUST', 'MUSE',
    'MASH', 'MASS', 'MASK', 'EASE', 'EAST', 'FAST', 'LAST', 'MAST', 'MEAT',
    'MEAD', 'MEAL', 'MEAN', 'MEAS', 'MESA'
])

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in WORDS:
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

# Find and print the path
path = find_path('WEAN', 'MESA')
if path:
    print(','.join(path))
else:
    print("No valid path found")