from collections import deque, defaultdict

# First, let's define a set of common 4-letter English words
words = set([
    'FOAL', 'FOOL', 'FOUL', 'FALL', 'FAIL', 'FOLD', 'FOOD', 'FOOT', 'FORD',
    'DROP', 'DRIP', 'PROP', 'CROP', 'DRAP', 'TRAP', 'TRIP', 'TROP', 'FROG',
    'COAL', 'COAT', 'GOAT', 'GOAL', 'ROAD', 'ROAM', 'FOAM', 'FROG', 'FROM',
    'FRAP', 'FLAP', 'FLIP', 'FLOP', 'FLOG', 'FROG', 'PROG', 'PROD'
])

def differ_by_one(word1, word2):
    return sum(1 for i in range(len(word1)) if word1[i] != word2[i]) == 1

def find_path(start, end):
    if start not in words or end not in words:
        return None
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        word, path = queue.popleft()
        if word == end:
            return path
            
        for next_word in words:
            if next_word not in visited and differ_by_one(word, next_word):
                visited.add(next_word)
                queue.append((next_word, path + [next_word]))
    
    return None

# Find the path
path = find_path('FOAL', 'DROP')
if path:
    print(','.join(path))
else:
    print("No path found")