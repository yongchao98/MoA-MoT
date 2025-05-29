from collections import deque, defaultdict

# Define a set of common 4-letter English words
WORDS = {
    'CELL', 'CELT', 'CENT', 'NECK', 'DECK', 'DELL', 'FELT', 'BELT', 'BENT', 'NEST',
    'NEXT', 'NEAT', 'NEAR', 'PELT', 'PENT', 'RENT', 'SENT', 'TENT', 'TELL', 'WELL',
    'WELT', 'WENT', 'BECK', 'PECK', 'PEEK', 'SEEK', 'WEEK', 'MELT', 'MELT', 'MEAT',
    'MEAD', 'MEAL', 'MEAN', 'LEAN', 'LEAD', 'LEND', 'LENS', 'LESS', 'MESS', 'MOSS',
    'MOST', 'MUST', 'MUSK', 'MASK', 'MARK', 'PARK', 'PACK', 'PICK', 'TICK', 'TACK',
    'TALK', 'TALL', 'TELL', 'TELL', 'SELL', 'SELF', 'SALT', 'HALT', 'HALL', 'HELL',
    'HELD', 'HEAD', 'HEAL', 'HEAR', 'HEIR', 'HELP', 'HELM', 'CELL', 'NELL', 'NERD',
    'NEED', 'FEED', 'FEEL', 'PEEL', 'PEER', 'PEAR', 'PEAK', 'PECK', 'NECK'
}

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
                new_path = list(path)
                new_path.append(neighbor)
                queue.append(new_path)
    
    return None

# Find the path and print the result
path = find_path('CELT', 'NECK')
if path:
    print(','.join(path))
else:
    print("No path found")