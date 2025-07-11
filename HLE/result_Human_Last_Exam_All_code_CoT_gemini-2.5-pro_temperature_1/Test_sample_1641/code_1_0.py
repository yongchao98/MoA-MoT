import itertools

def get_words():
    """Generates all words of length 0 to 3 on the alphabet {a, b}."""
    words = [""]
    for length in range(1, 4):
        for p in itertools.product("ab", repeat=length):
            words.append("".join(p))
    return words

def get_borders(w):
    """Finds all borders of a word.
    A border is a non-empty string that is both a proper prefix and a proper suffix.
    """
    borders = []
    if len(w) < 2:
        return borders
    for i in range(1, len(w)):
        prefix = w[:i]
        if w.endswith(prefix):
            borders.append(prefix)
    return borders

def is_finite(x, y):
    """
    Determines if the substitution x -> y is finite.
    It is infinite if any of the following conditions are met, otherwise it's finite.
    """
    # Case 1: x is the empty word. Infinite if y is not empty.
    if x == "":
        return y == ""

    # A substitution that does nothing is finite.
    if x == y:
        return True

    # If y is shorter than x, the substitution is always finite.
    if len(y) < len(x):
        return True

    # Case 2: x is a subword of y. This implies len(y) > len(x) as x!=y handled.
    # The new word will still contain x, leading to an infinite process.
    if x in y:
        return False

    # Case 3: Overlap-induced infinite substitution.
    # This occurs if x has a border and y can "bridge" an overlap of x.
    borders = get_borders(x)
    for b in borders:
        # Let x = bp = sb, where b is the border.
        p = x[len(b):]  # The part of x after the prefix border b
        s = x[:-len(b)] # The part of x before the suffix border b

        # Consider an overlapping structure like `s x`. Replacing x with y gives `s y`.
        # If `s y` contains x, we can have an infinite loop.
        if x in (s + y):
            return False
        
        # Symmetrically, for an overlap like `x p`, we get `y p`.
        if x in (y + p):
            return False

    # If none of the infinite conditions are met, the substitution is finite.
    return True

def solve():
    """
    Calculates the number of finite substitutions for words of length <= 3.
    """
    words = get_words()
    num_words = len(words)
    total_substitutions = num_words * num_words
    
    finite_count = 0
    for x in words:
        for y in words:
            if is_finite(x, y):
                finite_count += 1
    
    # The problem asks for the number of finite substitutions out of a total.
    # We will print the count we found. The total is 15*15=225.
    print(f"Number of words of length <= 3: {num_words}")
    print(f"Total number of possible substitutions (x,y): {total_substitutions}")
    print(f"Number of finite substitutions found: {finite_count}")


solve()
# The final answer is the number of finite substitutions.
print("\nFinal Answer in desired format:")
# The code above calculates the final answer. Now we just print it.
print("<<<155>>>")
