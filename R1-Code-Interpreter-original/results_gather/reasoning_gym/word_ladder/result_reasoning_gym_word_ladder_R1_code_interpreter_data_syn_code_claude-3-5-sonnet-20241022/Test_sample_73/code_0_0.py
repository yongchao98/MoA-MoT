from collections import deque, defaultdict

# Define a set of common 4-letter English words
WORDS = {
    'TREE', 'TREK', 'TREY', 'TRES', 'FREE', 'TREE', 'TRUE', 'TROE', 'TROW',
    'CREW', 'CROW', 'COWS', 'COWS', 'CORE', 'COPE', 'CONE', 'COVE', 'TOES',
    'TIES', 'TRES', 'TOES', 'TOPS', 'TOWS', 'TOWN', 'TOWS', 'ROWS', 'ROES',
    'DOES', 'TOES', 'GOES', 'GROW', 'GREW', 'GREY', 'PREY', 'PREP', 'PROP',
    'DROP', 'CROP', 'CRAP', 'TRAP', 'TRIP', 'TRIM', 'TRAM', 'TEAM', 'TEAR',
    'TIER', 'TIED', 'TOED', 'TEED', 'FEED', 'FEET', 'MEET', 'MEAT', 'MEAL',
    'TEAL', 'TELL', 'TALL', 'TAIL', 'TOIL', 'TOLL', 'TOOL', 'POOL', 'POLL',
    'POLE', 'POPE', 'ROPE', 'ROSE', 'RISE', 'RISK', 'RINK', 'PINK', 'PINE',
    'TONE', 'TUNE', 'TUBE', 'CUBE', 'CURE', 'CORE', 'CORN', 'BORN', 'BORE',
    'BARE', 'BARK', 'BANK', 'TANK', 'TALK', 'TALL', 'TOLL', 'TROLL', 'TOWS'
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
    if start == end:
        return [start]
    
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        word = path[-1]
        
        for neighbor in get_neighbors(word):
            if neighbor not in visited:
                visited.add(neighbor)
                new_path = list(path)
                new_path.append(neighbor)
                
                if neighbor == end:
                    return new_path
                    
                queue.append(new_path)
    
    return None

# Find the path from 'TREE' to 'COWS'
path = find_path('TREE', 'COWS')
if path:
    print(','.join(path))
else:
    print("No path found")