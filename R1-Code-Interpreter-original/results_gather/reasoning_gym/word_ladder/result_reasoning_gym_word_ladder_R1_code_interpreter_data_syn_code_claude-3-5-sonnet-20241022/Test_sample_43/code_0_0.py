from collections import deque, defaultdict

# Common 4-letter words
words = '''TUMS TOMB TOME TAME TAPE TAPS TOPS TOES TIES TINS TENS DENS DENT
TENT TERM TEAM TEAL TELL TALL TALE TAKE TAME SAME SOME SOLE SALE PALE PALM
PALL POLL POOL COOL COOK TOOK TOCK TICK TACK PACK PARK PERK PEEK PEEL PEER
BEER BEAR BEAT BELT BENT BEND BAND BOND BONE TONE TUNE TURN TERN TEEN SEEN
SEED SEND SAND HAND HARD HARM WARM WARD WORD WORK PORK PERK PEEK PEEK PEEL
FEEL FELL FILL PILL PILE PINE PING WING WINE WIRE WERE WEAR WEAK PEAK PEAR
TEAR TEAT THAT THAN THIN THIS THIS THUS THUD THEM THEN WHEN WHAT PREP'''

def build_graph(words):
    # Create graph where words are connected if they differ by one letter
    graph = defaultdict(list)
    word_set = set(words.split())
    
    for word in word_set:
        # Generate all possible one-letter variations
        for i in range(4):
            for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                new_word = word[:i] + c + word[i+1:]
                if new_word in word_set and new_word != word:
                    graph[word].append(new_word)
    return graph

def find_path(start, end, graph):
    queue = deque([[start]])
    visited = {start}
    
    while queue:
        path = queue.popleft()
        word = path[-1]
        
        if word == end:
            return path
            
        for neighbor in graph[word]:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(path + [neighbor])
    
    return None

# Find the path
graph = build_graph(words)
path = find_path('TUMS', 'PREP', graph)

if path:
    print(','.join(path))
else:
    print("No path found")