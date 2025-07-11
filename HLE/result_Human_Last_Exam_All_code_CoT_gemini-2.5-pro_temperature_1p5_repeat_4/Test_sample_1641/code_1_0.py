import itertools

def get_words():
    """Generates all words of length 0 to 3 on the alphabet {a, b}."""
    words = [""]
    for length in range(1, 4):
        for p in itertools.product(['a', 'b'], repeat=length):
            words.append("".join(p))
    return words

def get_proper_prefixes(word):
    """Returns a set of non-empty proper prefixes of a word."""
    return {word[:i] for i in range(1, len(word))}

def get_proper_suffixes(word):
    """Returns a set of non-empty proper suffixes of a word."""
    return {word[i:] for i in range(1, len(word))}

def is_finite(x, y):
    """
    Determines if the substitution x -> y is finite.
    """
    # Rule 1: Substitution from an empty string is infinite.
    if not x:
        return False

    # Rule 2: If y is shorter than x, the substitution is always finite.
    if len(y) < len(x):
        return True

    # Rule 3: If y contains x, it's infinite.
    if x in y:
        return False

    # Rule 4: Check for procreative overlaps.
    # Case a: A proper suffix of x is a proper prefix of y.
    proper_suffixes_x = get_proper_suffixes(x)
    proper_prefixes_y = get_proper_prefixes(y)
    if not proper_suffixes_x.isdisjoint(proper_prefixes_y):
        return False

    # Case b: A proper prefix of x is a proper suffix of y.
    proper_prefixes_x = get_proper_prefixes(x)
    proper_suffixes_y = get_proper_suffixes(y)
    if not proper_prefixes_x.isdisjoint(proper_suffixes_y):
        return False

    # If none of the infinite conditions are met, it's finite.
    return True

def solve():
    """
    Counts the number of finite substitutions for words of length <= 3.
    """
    words = get_words()
    
    counts = {
        "|y|<|x|": 0,
        "|x|=1, |y|>=1": 0,
        "|x|=2, |y|>=2": 0,
        "|x|=3, |y|=3": 0,
    }

    total_finite_count = 0
    for x in words:
        for y in words:
            if is_finite(x, y):
                total_finite_count += 1
                if len(y) < len(x):
                    counts["|y|<|x|"] += 1
                elif len(x) == 1 and len(y) >= 1:
                    counts["|x|=1, |y|>=1"] += 1
                elif len(x) == 2 and len(y) >= 2:
                     counts["|x|=2, |y|>=2"] += 1
                elif len(x) == 3 and len(y) == 3:
                     counts["|x|=3, |y|=3"] += 1
    
    print("The number of finite substitutions is the sum of counts from different categories:")
    
    parts = []
    for category, count in counts.items():
        if count > 0:
            parts.append(str(count))
            
    equation = " + ".join(parts)
    print(f"Total finite substitutions = {equation} = {total_finite_count}")


solve()
print("\n<<<88>>>")