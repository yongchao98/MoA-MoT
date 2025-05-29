from collections import deque, defaultdict

# Create a set of valid 4-letter English words
words = set('''TECH TACH TACK TECK DECK NECK TEAK TEAM TEAR TEAT TELL TEAL TEEN TEED TEES TENT TERM TERN TEST TEXT THEM THEN THEE TREE TREK TREY LEET LEFT LENT LEST LETS LEAD LEAF LEAK LEAN LEAP LEAS FEAT NEAT MEAT HEAT BEAT BELT BENT BEST LEST LESS LENS LEND LEAD'''.split())

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

# Find the path and print result
path = find_path('TECH', 'LEET')
if path:
    print(','.join(path))
else:
    print("No valid path found")