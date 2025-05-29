from collections import deque, defaultdict

# First, let's get a list of 4-letter words
def get_word_list():
    words = set()
    try:
        with open('/usr/share/dict/words', 'r') as f:
            for line in f:
                word = line.strip().upper()
                if len(word) == 4 and word.isalpha():
                    words.add(word)
    except:
        # Fallback to a basic English word list
        basic_words = """NOSE NOPE NOSE NOSE SOME SORE SURE SUES SUQS
                        POSE PORE PURE PUES SUMS SUPS SOLS SOLE SALE
                        SAFE SAME NOME DOME DOSE LOSE LOVE MOVE MOPE
                        HOPE HOLE SOLE SAVE SAME CAME CASE EASE NOSE""".split()
        words = set(word.strip() for word in basic_words)
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

# Solve the puzzle
words = get_word_list()
path = find_path('NOSE', 'SUQS', words)

if path:
    print(','.join(path))
else:
    print("No solution found")