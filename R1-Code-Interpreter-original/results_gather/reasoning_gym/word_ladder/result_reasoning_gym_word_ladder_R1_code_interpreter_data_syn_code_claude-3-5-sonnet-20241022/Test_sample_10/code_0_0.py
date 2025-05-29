from collections import deque, defaultdict

def get_word_list():
    # Create a set of 4-letter words
    words = set()
    try:
        with open('/usr/share/dict/words', 'r') as f:
            for line in f:
                word = line.strip().upper()
                if len(word) == 4 and word.isalpha():
                    words.add(word)
    except:
        # Fallback to a basic set of common 4-letter words
        words = {'HAND', 'HARD', 'HEAD', 'HEAL', 'HEAR', 'HEAT', 'HELD', 'HELP', 'HERD',
                'HERE', 'HERO', 'HOLD', 'HOLE', 'LEAD', 'LEAF', 'LEAN', 'LEAP', 'LEND',
                'LESS', 'LEST', 'LETS', 'FLED', 'FLEE', 'FLEW', 'FLEX', 'GLEN', 'GLEE'}
    return words

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            neighbors.append(new_word)
    return neighbors

def find_path(start, end, words):
    if start not in words or end not in words:
        return None
    
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        word = path[-1]
        
        if word == end:
            return path
            
        for neighbor in get_neighbors(word):
            if neighbor in words and neighbor not in visited:
                visited.add(neighbor)
                new_path = list(path)
                new_path.append(neighbor)
                queue.append(new_path)
    
    return None

# Find the path
words = get_word_list()
path = find_path('HAND', 'GLEE', words)

if path:
    print(','.join(path))
else:
    print("No valid path found")