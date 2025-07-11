def generate_words():
    """Generates all words of length 0 to 3 on the alphabet {a, b}."""
    words = ['']
    for length in range(1, 4):
        for i in range(2**length):
            # Create binary strings and map them to 'a' and 'b'
            bin_str = bin(i)[2:].zfill(length)
            words.append(bin_str.replace('0', 'a').replace('1', 'b'))
    return words

def get_borders(x):
    """Finds all non-empty words that are both a proper prefix and a suffix of x."""
    borders = []
    for i in range(1, len(x)):
        if x[:i] == x[-i:]:
            borders.append(x[:i])
    return borders

def is_finite(x, y):
    """
    Checks if the substitution x -> y is finite based on string rewriting theory.
    """
    # 1. If x is the empty string, the substitution is infinite.
    if not x:
        return False

    # 2. If y is shorter than x, the substitution is always finite.
    if len(y) < len(x):
        return True

    # 3. If y contains x as a subword, the substitution is infinite.
    if x in y:
        return False

    # 4. Check for infinite cases due to overlaps.
    borders = get_borders(x)
    for u in borders:
        # s is the prefix of x when u is the suffix: x = s + u
        s = x[:-len(u)]
        # If s followed by y creates a new x, it's infinite.
        if x in (s + y):
            return False
            
    # If none of the infinite conditions are met, it's finite.
    return True

def solve():
    """
    Counts the number of finite substitutions (x,y) for words of length <= 3.
    """
    words = generate_words()
    
    # Initialize counts for different categories
    finite_len_lt = 0
    finite_len_ge = 0
    
    # Check all pairs (x, y)
    for x in words:
        for y in words:
            if is_finite(x, y):
                if len(y) < len(x):
                    finite_len_lt += 1
                else:
                    finite_len_ge += 1

    total_finite = finite_len_lt + finite_len_ge
    
    print("Plan:")
    print("1. All 15 cases where x is the empty word are infinite.")
    print("2. Count finite cases where len(y) < len(x).")
    print("3. Count finite cases where len(y) >= len(x) and x is not empty.")
    print("\nCalculation:")
    print(f"Number of finite substitutions where len(y) < len(x): {finite_len_lt}")
    print(f"Number of finite substitutions where len(y) >= len(x): {finite_len_ge}")
    print("\nTotal finite substitutions:")
    print(f"{finite_len_lt} + {finite_len_ge} = {total_finite}")

solve()