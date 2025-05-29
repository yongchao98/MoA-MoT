# Define a set of common 4-letter words
words = {
    'LINE', 'LONE', 'LANE', 'LANG', 'LONG', 'LUNG', 'LOSE', 'LOGS', 'FOGS',
    'FINE', 'FIRE', 'FORE', 'FORK', 'FOLK', 'FOND', 'FEND', 'FIND', 'FING',
    'LING', 'LINK', 'LEND', 'LENS', 'LESS', 'LOSS', 'LOWS', 'LAWS', 'LAGS',
    'MINE', 'MINT', 'MIND', 'MEND', 'MENS', 'MONS', 'MOSS', 'MISS', 'MIST',
    'PINE', 'PING', 'PONG', 'POND', 'POLE', 'POPE', 'POPS', 'POGS',
    'WINE', 'WING', 'WONG', 'WONG', 'WOES', 'GOES', 'FOES'
}

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
    
    queue = [(start, [start])]
    visited = {start}
    
    while queue:
        word, path = queue.pop(0)
        
        if word == end:
            return path
            
        for neighbor in get_neighbors(word):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    
    return None

# Find path from LINE to FOGS
path = find_path('LINE', 'FOGS')
if path:
    print(','.join(path))
else:
    print("No path found")