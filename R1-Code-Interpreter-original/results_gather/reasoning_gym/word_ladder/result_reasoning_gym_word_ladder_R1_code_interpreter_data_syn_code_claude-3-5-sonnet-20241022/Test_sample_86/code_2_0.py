from collections import deque

# Carefully curated list focusing on possible connecting words
words = set('''
LAST CAST CASE CASH CASK COST COT COAT CHAT THAT THAN THEN THEE TREE TORE TONE
TINE TIME TAME TALE TALK WALK WALL WELL WERE WIRE WISE WIPE WINE WING WINK SINK
SINK SANK SAND BAND BEND BENT BEAT BOAT COAT COST LOST LOFT LIFT LIFE LIST LAST
CAST CASH CASE CARE CORE COPE CODE CODA SODA SOMA SOME COME COMB TOMB TONE TORE
FORE FORK WORK WORD WARD CARD CORD COLD HOLD HOLE HOPE HOSE LOSE LOVE MOVE MODE
MADE MALE MILE MINE MINT HINT HUNT HURT CURE PURE PORE POSE ROSE ROVE COVE COLE
COLT COAT COST CAST CASE CAVE COVE CODE CONE COPE CORE CORN BORN BORE BARE BASE
CASE CASH DASH DISH FISH FIST FAST FACT FACE FAME GAME GATE HATE HAVE HIVE HIRE
HIDE RIDE RISE ROSE POSE PORE PURE CURE CORE COKE CAKE LAKE LIKE LIVE LOVE MOVE
MOLE POLE POPE COPE CAPE CASE CASH CLASH FLASH FLESH FRESH FROST FRONT COUNT COURT
COAST CAST LAST LOST COST COSH COHO
'''.split())

def find_path(start, end, words):
    if start not in words or end not in words:
        return None
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        word, path = queue.popleft()
        if word == end:
            return path
            
        # Generate all possible words that differ by one letter
        for i in range(len(word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                new_word = word[:i] + c + word[i+1:]
                if new_word in words and new_word not in visited:
                    visited.add(new_word)
                    new_path = path + [new_word]
                    if new_word == end:
                        return new_path
                    queue.append((new_word, new_path))
    
    return None

# Find the path
path = find_path('LAST', 'COHO', words)
if path:
    print(','.join(path))
else:
    print("No path found")