from collections import deque

# Create a very specific set of 4-letter words focusing on possible paths
words = set('''
    BOON BOOK BOOT BOUT BOUT BOUT BOAT BEAT BRAT BRAT BART BARD BARE BORE BONE
    BOND BORN BARN BARK DARK DART PART PORT FORT FORM FROM ARTS ARTY
    BOON BORN BURN TURN TORN TORE TARE TART TARN YARN YARD CARD CART PART PERT PEAT PRAT BRAT
    BOON SOON SOOT SORT PORT PART ARTS ARTY
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
path = find_path('BOON', 'ARTY')
if path:
    print(','.join(path))
else:
    print("No path found")