import collections

def generate_words(max_len):
    """Generates all words of length up to max_len on alphabet {a, b}."""
    words = {""}
    for i in range(1, max_len + 1):
        new_words = set()
        for word in words:
            if len(word) == i - 1:
                new_words.add(word + 'a')
                new_words.add(word + 'b')
        words.update(new_words)
    return sorted(list(words), key=len)

def get_borders(x):
    """Finds all borders of a word x.
    A border is a non-empty word that is both a proper prefix and a proper suffix.
    """
    borders = set()
    for i in range(1, len(x)):
        prefix = x[:i]
        suffix = x[len(x)-i:]
        if prefix == suffix:
            borders.add(prefix)
    return borders

def is_finite(x, y):
    """Determines if the substitution x -> y is finite."""
    # Rule 1: x is empty.
    if not x:
        return False
    
    # Rule 2: len(y) < len(x).
    if len(y) < len(x):
        return True
    
    # Rule 3: len(y) == len(x).
    if len(y) == len(x):
        return x != y
        
    # Rule 4: len(y) > len(x).
    if x in y:
        return False
        
    borders = get_borders(x)
    if not borders:
        # If there are no borders, no overlap can create x.
        return True
    
    for b in borders:
        # x = s + b = b + p. Note: for python slicing, if x = b + p, then p = x[len(b):]
        s = x[:-len(b)]
        p = x[len(b):]
        
        # Check for overlap resubstitution
        if x in (s + y) or x in (y + p):
            return False
            
    return True

def solve():
    """
    Counts the number of finite substitutions for words of length <= 3.
    The prompt mentions 255 substitutions, which suggests 16*16 - 1 pairs.
    However, words of length <= 3 on {a,b} number 1+2+4+8 = 15.
    This leads to 15*15 = 225 possible pairs. This code solves for these 225 pairs.
    """
    words = generate_words(3)
    
    total_pairs = len(words) * len(words)
    finite_count = 0
    
    # To store counts for the final equation
    counts = collections.defaultdict(int)

    for x in words:
        for y in words:
            if is_finite(x, y):
                finite_count += 1
                if len(y) < len(x):
                    counts['len_less'] += 1
                elif len(y) == len(x):
                    counts['len_equal'] += 1
                else: # len(y) > len(x)
                    counts['len_greater'] += 1

    c_less = counts['len_less']
    c_equal = counts['len_equal']
    c_greater = counts['len_greater']
    
    print(f"Number of finite substitutions where len(y) < len(x): {c_less}")
    print(f"Number of finite substitutions where len(y) == len(x): {c_equal}")
    print(f"Number of finite substitutions where len(y) > len(x): {c_greater}")
    print(f"Total finite substitutions: {c_less} + {c_equal} + {c_greater} = {finite_count}")


solve()
<<<156>>>