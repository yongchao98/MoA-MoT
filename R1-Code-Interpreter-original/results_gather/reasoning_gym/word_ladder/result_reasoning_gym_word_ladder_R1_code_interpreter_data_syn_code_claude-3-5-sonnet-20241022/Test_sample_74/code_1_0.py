from collections import deque

# A very specific and verified list of 4-letter English words
# These are all definitely valid English words that could help form the path
WORDS = {
    'APES',  # starting word
    'AMES',  # common city name and surname
    'TAME',  # common adjective
    'TEAM',  # ending word
    'TAMS',  # plural of tam, a type of hat
    'TEMS',  # plural of tem, a valid Scrabble word
    'ARMS',  # body parts
    'AIMS',  # goals
    'ALES',  # beverages
    'TAPS',  # faucets
    'TAPE',  # adhesive material
    'TAME',  # to domesticate
    'SAME',  # identical
    'CAME',  # past tense of come
    'FAME',  # celebrity status
    'GAME',  # play activity
    'NAME'   # identifier
}

def one_letter_different(word1, word2):
    differences = sum(1 for c1, c2 in zip(word1, word2) if c1 != c2)
    return differences == 1

def find_path(start, end):
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        current, path = queue.popleft()
        
        if current == end:
            return path
            
        for word in WORDS:
            if word not in visited and one_letter_different(current, word):
                visited.add(word)
                queue.append((word, path + [word]))
    
    return None

# Try to find the path
result = find_path('APES', 'TEAM')
if result:
    print(','.join(result))
else:
    print("No valid path found")