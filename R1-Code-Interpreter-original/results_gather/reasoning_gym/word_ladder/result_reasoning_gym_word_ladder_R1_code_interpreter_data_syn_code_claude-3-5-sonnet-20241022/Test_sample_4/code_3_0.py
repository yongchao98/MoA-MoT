from collections import deque

# Create a focused set of 4-letter words that are more likely to form a path
words = set('''
    BOON BOOK BOOT BOOM ROOM ROAM ROAD READ REAL REAR BEAR BEER BEEF REEF REEL
    REAL DEAL DEAD DEAR DEER DEEP BEEP BEEN SEEN SEED SEND SAND BAND BANK BARK
    BARN BORN BORE BARE CARE CART CARD CARS EARS EARN TURN TORN TOWN DOWN DOOM
    ROOM ROAM FOAM FEAR HEAR HEAT NEAT NEXT TEXT TENT WENT WANT WART WARD WARM
    WORM WORN BURN TURN TERN TERM TEAM BEAM BEAR HEAR HEAT HEAD HARD CARD CARE
    DARE DARK PARK PART ARTS ARTY
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