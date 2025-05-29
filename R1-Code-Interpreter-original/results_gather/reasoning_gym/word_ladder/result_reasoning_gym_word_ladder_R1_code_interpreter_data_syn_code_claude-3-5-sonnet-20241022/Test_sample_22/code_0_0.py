from collections import deque, defaultdict

# First, let's create a set of valid 4-letter words
def get_words():
    words = set()
    try:
        with open('/usr/share/dict/words', 'r') as f:
            for line in f:
                word = line.strip().upper()
                if len(word) == 4 and word.isalpha():
                    words.add(word)
    except:
        # Fallback to a basic set of common 4-letter words
        words = {'PSIS', 'PEAS', 'TEAS', 'TEAT', 'PEST', 'PAST', 'PASS', 'TESS', 'TEST', 'TEAM', 'TEAR', 'PEAT'}
    return words

def get_neighbors(word):
    neighbors = []
    for i in range(4):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            neighbors.append(new_word)
    return neighbors

def find_path(start, end, words):
    if start not in words or end not in words:
        return None
    
    queue = deque([[start]])
    seen = {start}
    
    while queue:
        path = queue.popleft()
        word = path[-1]
        
        if word == end:
            return path
            
        for neighbor in get_neighbors(word):
            if neighbor in words and neighbor not in seen:
                queue.append(path + [neighbor])
                seen.add(neighbor)
    
    return None

# Solve the word ladder
words = get_words()
path = find_path('PSIS', 'TEAT', words)

if path:
    print(','.join(path))
else:
    print("No solution found")