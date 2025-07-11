import itertools

def solve():
    """
    Calculates the number of finite substitutions (x,y) for words of length <= 3
    on the alphabet {a,b}.
    """
    max_len = 3
    words = [""]
    alphabet = ['a', 'b']
    for length in range(1, max_len + 1):
        for p in itertools.product(alphabet, repeat=length):
            words.append("".join(p))

    def is_infinite_growing(x, y):
        """
        Checks for infinite substitution when len(y) > len(x).
        This happens if a new 'x' can be created.
        """
        # Case 1: x is a substring of y
        if x in y:
            return True
        
        # Case 2: x is formed by context. For every split x=uv,
        # check if x is a substring of yu or vy.
        for i in range(1, len(x)):
            u = x[:i]  # prefix
            v = x[i:]  # suffix
            # A new x is created if it's found in y concatenated with
            # a prefix of x, or a suffix of x concatenated with y.
            if x in (y + u) or x in (v + y):
                return True
        return False

    finite_count = 0
    total_substitutions = len(words) * len(words)

    for x in words:
        for y in words:
            is_finite = False
            
            # Case: x is the empty word
            if not x:
                if not y:
                    is_finite = True
            # Case: substitution shortens the string
            elif len(y) < len(x):
                is_finite = True
            # Case: substitution preserves length
            elif len(y) == len(x):
                if x != y:
                    is_finite = True
            # Case: substitution grows the string
            else: # len(y) > len(x)
                if not is_infinite_growing(x, y):
                    is_finite = True
            
            if is_finite:
                finite_count += 1
                
    print(f"Total number of substitutions considered: {total_substitutions}")
    print(f"Number of finite substitutions: {finite_count}")

solve()
print("<<<158>>>")