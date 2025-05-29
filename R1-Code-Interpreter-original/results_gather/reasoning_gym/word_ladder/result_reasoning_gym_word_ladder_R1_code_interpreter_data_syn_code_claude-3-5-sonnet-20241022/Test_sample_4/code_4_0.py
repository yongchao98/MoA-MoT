from collections import deque

# Create a comprehensive set of 4-letter words
words = set('''
    BOON BOOK BOOM BOOT BOUT BOLT BONE BORE BORN BURN BARN BARK BARE CARE CART
    CAST COST COAT BOAT BEAT BENT BEST REST PEST PART PANT PAST CAST CART CARE
    DARE DARK PARK PART PANT PINT PINE PILE POLE PORE PURE PARE PART ARTS ARTY
    BORN WORN WARN WART WARD CARD CART CAST COST COT COT CAT CART PART PARE
    PARK DARK DART DIRT DIRE FIRE FARE CARE CART CAST COST COAT BOAT BEAT BENT
    BERT PERT PART ARTS ARTY
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