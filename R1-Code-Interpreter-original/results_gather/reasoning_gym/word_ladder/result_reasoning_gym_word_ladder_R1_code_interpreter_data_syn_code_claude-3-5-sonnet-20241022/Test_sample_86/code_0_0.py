from collections import deque, defaultdict

# First, let's create a set of common 4-letter words
words = set('''
LAST LOST COST CAST FAST PAST MAST MOST HOST POST COST COSY COSH CASH CASE
LASH LACK LACE LATE LAKE LAME DAME DARE CARE CORE COPE CODE CODE COME COMB
COAT COST COLT COLD HOLD HOLE HOME SOME SAME SAVE HAVE HOVE HOSE HOLE HOLY
COHO
'''.split())

def diff_by_one(word1, word2):
    return sum(1 for i in range(len(word1)) if word1[i] != word2[i]) == 1

def find_path(start, end, words):
    if start not in words or end not in words:
        return None
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        word, path = queue.popleft()
        if word == end:
            return path
            
        # Try changing each position
        for i in range(len(word)):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                new_word = word[:i] + c + word[i+1:]
                if new_word in words and new_word not in visited:
                    visited.add(new_word)
                    queue.append((new_word, path + [new_word]))
    
    return None

# Find the path
path = find_path('LAST', 'COHO', words)
if path:
    print(','.join(path))
else:
    print("No path found")