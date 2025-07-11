import collections

def is_finite_substitution(x, y):
    """
    Determines if the substitution x -> y is finite based on established criteria.
    """
    # Case 1: x is the empty word
    if not x:
        return not y

    # Case 2: Length-decreasing substitution
    if len(y) < len(x):
        return True

    # Case 3: Length-preserving substitution
    if len(y) == len(x):
        # A rule is infinite if x and y are conjugates (and x != y).
        # x and y are conjugates if x = uv and y = vu for non-empty u,v.
        # This is equivalent to checking if y is a non-trivial cyclic shift of x.
        if x == y:
            return True
        if len(x) > 0 and y in (x + x):
            return False # They are conjugates, so infinite.
        return True # Not conjugates, so finite.

    # Case 4: Length-increasing substitution
    if len(y) > len(x):
        # Infinite if y contains x as a subword.
        if x in y:
            return False
        
        # Infinite if an overlap can create a new x.
        # This happens if a prefix of x followed by y, or y followed by a suffix of x, contains x.
        n = len(x)
        # Check `prefix + y` for all proper non-empty prefixes of x
        for i in range(1, n):
            prefix = x[:i]
            if x in (prefix + y):
                return False
        
        # Check `y + suffix` for all proper non-empty suffixes of x
        for i in range(1, n):
            suffix = x[i:]
            if x in (y + suffix):
                return False

        # If no conditions for infinitude are met, it's finite.
        return True
        
    return False # Should not be reached

def solve():
    """
    Generates all words, iterates through all pairs (x,y), counts the finite ones, and prints the result.
    """
    max_len = 3
    alphabet = ['a', 'b']
    
    # Generate all words of length 0 to max_len
    words = [""]
    queue = collections.deque([""])
    while queue:
        current_word = queue.popleft()
        if len(current_word) < max_len:
            for char in alphabet:
                new_word = current_word + char
                words.append(new_word)
                queue.append(new_word)

    finite_count = 0
    for x in words:
        for y in words:
            if is_finite_substitution(x, y):
                finite_count += 1
                
    print("The number of finite substitutions is:")
    print(finite_count)

solve()