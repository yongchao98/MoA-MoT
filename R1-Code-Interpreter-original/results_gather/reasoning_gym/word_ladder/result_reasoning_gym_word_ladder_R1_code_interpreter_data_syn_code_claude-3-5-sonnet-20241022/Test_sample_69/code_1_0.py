from collections import deque

# Expanded word list including more common 4-letter words
words = '''MESS MUSS MUSE SUES SEES SEED MEET MUST MIST MISS
          MASS MOSS LESS LOSE LOST LAST FAST PAST PASS BASS
          BASE CASE CAST COST COAT BOAT BEAT MEAT MEAL MEAD
          LEAD LOAD ROAD READ REAL SEAL SEAT SENT TENT TENS
          TEEN SEEN SEEM STEM STEP STOP SHOP SHIP SLIP SLIM
          SLAM SLUM GLUM GLUE BLUE BLUR SUER SUED SESS MEWS
          MUSS FUSS FESS TESS TOSS MOSS MOPS MAPS MATS MATS
          SETS SITS SINS PINS PENS PETS PUTS PUSS PULS PALS'''.split()

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in word_set and new_word != word:
                neighbors.append(new_word)
    return neighbors

# Create word set for O(1) lookup
word_set = set(words)

def find_ladder(start, end):
    if start not in word_set or end not in word_set:
        return None
    
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        current = path[-1]
        
        if current == end:
            return path
        
        for i in range(len(current)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                next_word = current[:i] + c + current[i+1:]
                if next_word in word_set and next_word not in visited:
                    visited.add(next_word)
                    new_path = list(path)
                    new_path.append(next_word)
                    queue.append(new_path)
    return None

# Find and print the path
result = find_ladder('SUES', 'MESS')
if result:
    print(','.join(result))
else:
    print("No solution found")