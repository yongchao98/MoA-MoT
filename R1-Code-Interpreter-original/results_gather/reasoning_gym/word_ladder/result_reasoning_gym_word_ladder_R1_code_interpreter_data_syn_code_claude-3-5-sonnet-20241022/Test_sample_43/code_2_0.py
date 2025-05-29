def can_transform(word1, word2):
    return sum(1 for c1, c2 in zip(word1, word2) if c1 != c2) == 1

# Very carefully curated list ensuring connectivity
words = '''
TUMS TUBS TUTS PUTS PUNS PENS PERS PEPS PREP
TUMS TUNS TENS TERS TEPS TEMP TUMP PUMP PUMP
TUMS TUTS TOTS TOPS TOPS TAPS TIPS TRIPS PREP
'''.strip().split()

def find_all_paths(start, end, word_list, max_length=10):
    def dfs(current_word, path, visited):
        if len(path) > max_length:
            return []
        if current_word == end:
            return [path]
        
        paths = []
        for word in word_list:
            if word not in visited and can_transform(current_word, word):
                new_paths = dfs(word, path + [word], visited | {word})
                paths.extend(new_paths)
        return paths

    all_paths = dfs(start, [start], {start})
    return sorted(all_paths, key=len)

# Find all possible paths
paths = find_all_paths('TUMS', 'PREP', words)

# Print the shortest path if found
if paths:
    print(','.join(paths[0]))
else:
    print("No path found")