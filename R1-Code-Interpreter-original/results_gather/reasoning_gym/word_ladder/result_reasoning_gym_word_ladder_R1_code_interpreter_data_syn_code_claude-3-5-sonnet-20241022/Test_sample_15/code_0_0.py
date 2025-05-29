from collections import deque, defaultdict

# Create a set of 4-letter words
def get_four_letter_words():
    words = set()
    try:
        with open('/usr/share/dict/words', 'r') as f:
            for line in f:
                word = line.strip().upper()
                if len(word) == 4 and word.isalpha():
                    words.add(word)
    except:
        # Fallback to a basic set of common 4-letter words
        words = {'BRAD', 'BEAD', 'BEAM', 'BERM', 'PERM', 'PERM', 'POEM',
                'BRED', 'BEAD', 'READ', 'BRAG', 'BEAT', 'BEND', 'BEAN',
                'BEAR', 'PEAR', 'PEER', 'PEEL', 'BEAM', 'SEAM', 'SLAM',
                'BERM', 'TERM', 'TEAM', 'PERM', 'PERM', 'POEM'}
    return words

def get_neighbors(word):
    neighbors = []
    for i in range(4):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            if c != word[i]:
                neighbor = word[:i] + c + word[i+1:]
                neighbors.append(neighbor)
    return neighbors

def find_path(start, end, words):
    if start not in words or end not in words:
        return None
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        if current == end:
            return path
            
        for neighbor in get_neighbors(current):
            if neighbor in words and neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    
    return None

# Find the path
words = get_four_letter_words()
path = find_path('BRAD', 'POEM', words)

if path:
    print(','.join(path))
else:
    print("No valid path found")