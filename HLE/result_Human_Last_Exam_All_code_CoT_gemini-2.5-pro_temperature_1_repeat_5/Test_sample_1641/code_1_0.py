import itertools

def is_finite(x, y):
    """
    Determines if the substitution x -> y is finite.
    """
    # Case 0: x is the empty word. This is always infinite.
    if not x:
        return False

    # Case 1: Shrinking substitution. Always finite.
    if len(y) < len(x):
        return True

    # Case 2: Length-preserving substitution.
    if len(y) == len(x):
        # Infinite only if x and y are identical.
        return x != y

    # Case 3: Expanding substitution (len(y) > len(x)).
    # 3a: If y contains x, it's a simple infinite loop.
    if x in y:
        return False

    # 3b: Check for overlaps that allow infinite pumping.
    # This happens if a prefix of y is a suffix of x, or a suffix of y is a prefix of x.

    # Check if any non-empty prefix of x is a suffix of y.
    for i in range(1, len(x) + 1):
        if y.endswith(x[:i]):
            return False

    # Check if any non-empty suffix of x is a prefix of y.
    for i in range(1, len(x) + 1):
        if y.startswith(x[-i:]):
            return False
            
    # If none of the infinite conditions are met, it's finite.
    return True

def solve():
    """
    Solves the problem by generating words, checking each pair, and counting.
    """
    words = [""]
    alphabet = ['a', 'b']
    for length in range(1, 3 + 1):
        for p in itertools.product(alphabet, repeat=length):
            words.append("".join(p))
    
    total_substitutions = len(words) * len(words)
    finite_count = 0
    
    for x in words:
        for y in words:
            if is_finite(x, y):
                finite_count += 1
                
    print(f"Total number of words of length <= 3 is {len(words)}.")
    print(f"Total substitutions considered (15*15) is {total_substitutions}.")
    print(f"Out of these, the number of finite substitutions is {finite_count}.")

solve()
<<<150>>>