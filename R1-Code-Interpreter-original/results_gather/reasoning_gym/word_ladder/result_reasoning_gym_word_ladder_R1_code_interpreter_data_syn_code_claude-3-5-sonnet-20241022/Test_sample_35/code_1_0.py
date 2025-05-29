from collections import deque

# Very carefully curated dictionary focusing on common 4-letter words
WORDS = {
    # Core words we know we need
    'CUSP', 'OVER',
    # Common 4-letter words that might help form the bridge
    'COST', 'CAST', 'CASE', 'CAVE', 'COVE',
    'COPE', 'CORE', 'CARE', 'CURE', 'CUPS',
    'CUTS', 'CUSS', 'OUST', 'DUST', 'MUST',
    'MUSE', 'MAST', 'MOST', 'MOVE', 'MORE',
    'CURE', 'PURE', 'PORE', 'POSE', 'ROSE',
    'ROVE', 'CASK', 'CASH', 'CAST', 'COAT',
    'CULT', 'CURT', 'CART', 'PART', 'PORT',
    'FORT', 'FORE', 'FIRE', 'FIVE', 'GIVE',
    # Additional bridging words
    'CASP', 'COSP', 'CROP', 'CRAP', 'TRAP',
    'TRIP', 'TRIM', 'TRAM', 'TEAM', 'TEAR',
    'TIER', 'TIER', 'TOES', 'GOES', 'GORE',
    'SORE', 'SURE', 'CURE', 'CORE', 'COPE',
    # More possible transitions
    'CUST', 'GUST', 'OUST', 'OAST', 'EAST',
    'EASE', 'EVER', 'EVEN', 'OVEN', 'OPEN'
}

def hamming_distance(word1, word2):
    return sum(c1 != c2 for c1, c2 in zip(word1, word2))

def get_neighbors(word):
    return [w for w in WORDS if hamming_distance(word, w) == 1]

def find_path(start, end):
    if start not in WORDS or end not in WORDS:
        return None
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        
        if current == end:
            return path
            
        for neighbor in get_neighbors(current):
            if neighbor not in visited:
                visited.add(neighbor)
                new_path = path + [neighbor]
                queue.append((neighbor, new_path))
    
    return None

# Try multiple starting points
paths = []
result = find_path('CUSP', 'OVER')
if result:
    print(','.join(result))
else:
    # Try some alternative paths through known intermediate words
    intermediate_words = ['CUSS', 'OUST', 'COST', 'COVE']
    for mid in intermediate_words:
        path1 = find_path('CUSP', mid)
        if path1:
            path2 = find_path(mid, 'OVER')
            if path2:
                full_path = path1 + path2[1:]
                paths.append(full_path)
                print(','.join(full_path))
                break