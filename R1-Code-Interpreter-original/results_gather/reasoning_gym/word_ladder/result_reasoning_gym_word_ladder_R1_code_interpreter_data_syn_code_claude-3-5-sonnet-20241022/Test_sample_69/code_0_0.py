from collections import deque, defaultdict

# Common 4-letter English words
words = '''MESS MUSS MUSE SUES SEES SEED MEET MUST MIST MISS
          MASS MOSS LESS LOSE LOST LAST FAST PAST PASS BASS
          BASE CASE CAST COST COAT BOAT BEAT MEAT MEAL MEAD
          LEAD LOAD ROAD READ READ REAL SEAL SEAT SENT TENT
          TENS TEEN SEEN SEEM STEM STEP STOP SHOP SHIP SLIP
          SLIM SLAM SLUM GLUM GLUE BLUE BLUR SUER SUED'''.split()

# Create set for O(1) lookup
word_set = set(words)

def get_neighbors(word):
    neighbors = []
    for i in range(4):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in word_set and new_word != word:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    if start not in word_set or end not in word_set:
        return None
    
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

# Find and print the path
path = find_path('SUES', 'MESS')
if path:
    print(','.join(path))
else:
    print("No solution found")