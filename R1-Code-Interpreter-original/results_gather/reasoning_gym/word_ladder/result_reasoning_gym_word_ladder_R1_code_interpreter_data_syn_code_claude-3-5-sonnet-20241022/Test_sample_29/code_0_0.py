from collections import deque, defaultdict

# First, let's create a set of common 4-letter words
words = set('''
FEEL FEED FEET FELL HEEL PEEL REEL TEEL TEEN TEEM TEAM
TELL TEAL TEAR TERM THEM THEN MEET MELT MEED DEED DEEP
TEED REED REEF REEK PEEK PEER PEAR FEAR FEAT FEAT TEAT
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
            if neighbor not in visited:
                visited.add(neighbor)
                new_path = list(path)
                new_path.append(neighbor)
                queue.append(new_path)
    
    return None

# Find the path from FEEL to TEEM
path = find_path('FEEL', 'TEEM')
if path:
    print(','.join(path))
else:
    print("No path found")