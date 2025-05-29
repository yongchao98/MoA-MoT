from collections import deque

# Comprehensive list of 4-letter English words focusing on possible paths
words = '''CLAN CRAN CRAY CRAW CRAM CRAP CLAP CHAP CHAT THAT THAN THIN
          CLAN CLAD CHAD CHAR CHAP CHIP SHIP SLIP SLOP STOP STEP STEM
          CLAN CLAW FLAW FLAN PLAN PLUS PLUM PLUMP PUMP RUMP RUMP RUSH
          CLAN CRAN GRAN GRIN GAIN RAIN RUIN RUIN RUSH
          CLAN CRAN CRAP TRAP TRIP TRIM GRIM GRIN GAIN RAIN RUIN RUSH
          CLAN CLEAN CLEAR CLEAT CHEAT CHEST CHASE CHOSE CLOSE CLONE CLOWN
          CLAN CRAN GRAM GRIM GRIP DRIP DROP CROP CRAP TRAP TRAM TEAM
          CLAN PLAN PLAY PRAY GRAY GRAB CRAB CRAP TRAP TRIP DRIP GRIP
          CLAN CLAW CLAD GLAD GOAD GOLD GOLF GULF GULL GULL GUSH RUSH
          CLAN CRAN GRAN GRIN GRIT GRID GRIM GRUM GRUM GUSH RUSH'''.split()

def get_neighbors(word):
    neighbors = []
    for i in range(4):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in word_set and new_word != word:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        word, path = queue.popleft()
        for neighbor in get_neighbors(word):
            if neighbor == end:
                return path + [end]
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    return None

# Convert words list to set for O(1) lookup
word_set = set(words)

# Find the path
path = find_path('CLAN', 'RUSH')

if path:
    print(','.join(path))
else:
    print("No path found")