from collections import deque

# A more comprehensive list of 4-letter words that are more likely to be valid English words
words = set('''
REIN RAIN RAIN RAID RAIL RAIS REIS RIES RISE RISK RICK RICH RICE RACE RAGE RARE CARE CASE CASH CAST COST COAT CHAT THAT THIS THIN CHIN CHIS
REIN VEIN VAIN RAIN RAIS REIS RIES RISE WISE WIRE WERE WEAR WEAK WEEK SEEK SEEN SHIN THIN CHIN CHIS
REIN RAIN MAIN MAIL MALL MALE MAKE MARE MORE MOVE LOVE LIVE LIKE PIKE PINE PINE FINE FIRE FIVE GIVE GIVE GIBE GIBE THIS CHIS
REIN RAIS RAID SAID SAND SEND SEED SEEK PEEK PERK PORK FORK FORM FROM FROG FLOG FLAG FLAN PLAN PLAY PRAY GRAY GRAS GRIS CRIS CHIS
'''.split())

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in words and sum(1 for a, b in zip(word, new_word) if a != b) == 1:
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
                queue.append(path + [neighbor])
    return None

# Try to find the path
path = find_path('REIN', 'CHIS')
if path:
    print(','.join(path))
else:
    print("No valid path found")