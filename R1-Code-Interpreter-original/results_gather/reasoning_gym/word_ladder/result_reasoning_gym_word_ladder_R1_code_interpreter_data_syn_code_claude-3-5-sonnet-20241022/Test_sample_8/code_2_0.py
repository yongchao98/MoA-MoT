# Using a more comprehensive word list focusing on common 4-letter words
WORDS = {
    'TOWS', 'SOWS', 'SOLS', 'SOIL', 'SOUL', 'SOUR', 'SOAR', 'SOAK',
    'TOSS', 'TOES', 'TOPS', 'TONS', 'TINS', 'TIES', 'TOWS', 'TOWN',
    'SOWS', 'SOLS', 'SOIL', 'TOIL', 'TAIL', 'TOLL', 'TOOL', 'POOL',
    'SOOT', 'SHOT', 'SPOT', 'SPIT', 'SPIN', 'SPAN', 'SOAP', 'SOAR',
    'SOUR', 'TOUR', 'TOLL', 'TELL', 'TILL', 'TOIL', 'SOIL', 'SOUL',
    'SODS', 'SONS', 'SINS', 'SITS', 'SIPS', 'SOPS', 'TOPS', 'TAPS',
    'TOWS', 'SOWS', 'SODS', 'SOWS', 'SOLS', 'SOIL'
}

def get_one_letter_changes(word):
    result = []
    for i in range(4):  # for each position
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            if c != word[i]:
                new_word = word[:i] + c + word[i+1:]
                if new_word in WORDS:
                    result.append(new_word)
    return result

def find_word_ladder(start, end, max_depth=6):
    if start not in WORDS or end not in WORDS:
        return None
    
    visited = {start: [start]}
    current_words = {start}
    depth = 0
    
    while current_words and depth < max_depth:
        depth += 1
        next_words = set()
        
        for word in current_words:
            for new_word in get_one_letter_changes(word):
                if new_word not in visited:
                    visited[new_word] = visited[word] + [new_word]
                    if new_word == end:
                        return visited[new_word]
                    next_words.add(new_word)
        
        current_words = next_words
    
    return None

# Try to find the path
result = find_word_ladder('TOWS', 'SOIL')
if result:
    print(','.join(result))
else:
    print("No path found")