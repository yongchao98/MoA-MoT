from collections import deque, defaultdict

# Define a set of common 4-letter English words
WORDS = set([
    'FORK', 'BELL', 'FOLK', 'FALL', 'FELL', 'FILL', 'FULL', 'FAIL', 'FOIL',
    'FOLL', 'POLL', 'PELL', 'PILL', 'PULL', 'PALL', 'BALL', 'BILL', 'BULL',
    'BELT', 'BALK', 'BULK', 'BOLL', 'BOLL', 'BOLL', 'BOLL', 'BOLL', 'BOLL',
    'BORK', 'BARK', 'BANK', 'BACK', 'PACK', 'PICK', 'PORK', 'PERK', 'PEAK',
    'BEAK', 'BEAM', 'BEAR', 'BEAT', 'BELT', 'SELL', 'TELL', 'WELL', 'WALL',
    'WALK', 'TALK', 'TALL', 'TELL', 'CELL', 'CALL', 'HALL', 'HELL', 'HEAL',
    'FOLD', 'BOLD', 'BELD', 'HELD', 'HELP', 'FEEL', 'FEED', 'FEET', 'FORT',
    'FORM', 'FOAM', 'FEAR', 'BEAR', 'BEER', 'BEEN', 'BELL'
])

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in WORDS and new_word != word:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    if start not in WORDS or end not in WORDS:
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
                queue.append(path + [neighbor])
    
    return None

# Find the path and print result
path = find_path('FORK', 'BELL')
if path:
    print(','.join(path))
else:
    print("No valid path found")