from collections import deque, defaultdict

# Create a set of common 4-letter words
words = set('''
JIGS JAGS JOGS JETS JEST JESS JETS JIGS JOBS JOTS
JEST BEST TEST TEST TESS LESS LOSS TOSS BOSS BASS
JIGS JAGS BAGS BATS BITS SITS SETS JETS JEST
'''.split())

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in words:
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
                new_path = list(path)
                new_path.append(neighbor)
                queue.append(new_path)
    
    return None

# Find the path
path = find_path('JIGS', 'JEST')
if path:
    print(','.join(path))
else:
    print("No solution found")