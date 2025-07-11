import collections

def generate_words(alphabet, max_len):
    """Generates all words on the alphabet up to a maximum length."""
    words = [""]
    for i in range(1, max_len + 1):
        new_words = []
        for word in words:
            if len(word) == i - 1:
                for char in alphabet:
                    new_words.append(word + char)
        words.extend(new_words)
    # The above logic is a bit complex, let's simplify.
    words = [""]
    q = collections.deque([""])
    while q:
        curr = q.popleft()
        if len(curr) < max_len:
            for char in alphabet:
                next_word = curr + char
                words.append(next_word)
                q.append(next_word)
    return words

def is_infinite_overlap(x, y):
    """
    Checks for non-terminating overlaps for the case |y| > |x| > 0.
    A rule x -> y is infinite if:
    1. x is a subword of y.
    2. x = uv and y = svu for non-empty u, v, s.
    """
    if x in y:
        return True
    
    # Check for pumpable overlap: x=uv, y=svu
    for i in range(1, len(x)):
        u = x[:i]
        v = x[i:]
        # u and v must be non-empty, which is guaranteed by the loop range.
        
        # Check if y can be written as s + v + u
        # This means y must start with v and end with u.
        if y.startswith(v) and y.endswith(u):
            # y = v...u. The part in the middle is s.
            # We need len(y) > len(v) + len(u) = len(x) for s to be non-empty.
            if len(y) > len(x):
                return True
    return False

def solve():
    """
    Counts the number of finite substitutions for words of length <= 3.
    """
    alphabet = ['a', 'b']
    max_len = 3
    words = generate_words(alphabet, max_len)
    
    total_pairs = len(words) * len(words)
    finite_count = 0
    
    # Counters for each category of finite substitutions
    finite_len_y_lt_x = 0
    finite_len_y_eq_x = 0
    finite_len_y_gt_x = 0

    for x in words:
        for y in words:
            is_finite = False
            # Case: x is the empty word. Always infinite.
            if not x:
                is_finite = False
            # Case: |y| < |x|. Always finite.
            elif len(y) < len(x):
                is_finite = True
                finite_len_y_lt_x += 1
            # Case: |y| == |x|.
            elif len(y) == len(x):
                # Finite if x != y, infinite if x == y.
                if x != y:
                    is_finite = True
                    finite_len_y_eq_x += 1
                else: # x == y
                    is_finite = False
            # Case: |y| > |x|.
            else: # len(y) > len(x)
                # This is the complex case. It's infinite if there's a non-terminating overlap.
                if is_infinite_overlap(x, y):
                    is_finite = False
                else:
                    is_finite = True
                    finite_len_y_gt_x += 1
            
            if is_finite:
                finite_count += 1

    print(f"Number of finite substitutions when |y| < |x|: {finite_len_y_lt_x}")
    print(f"Number of finite substitutions when |y| = |x| (and x!=y): {finite_len_y_eq_x}")
    print(f"Number of finite substitutions when |y| > |x|: {finite_len_y_gt_x}")
    print(f"Total finite substitutions: {finite_len_y_lt_x} + {finite_len_y_eq_x} + {finite_len_y_gt_x} = {finite_count}")
    print(f"Out of {total_pairs} total substitutions, {finite_count} are finite.")
    print(f"<<<{finite_count}>>>")

solve()