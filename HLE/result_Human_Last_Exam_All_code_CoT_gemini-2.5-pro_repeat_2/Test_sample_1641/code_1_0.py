import collections

def is_conjugate(x, y):
    """Checks if x and y are non-equal conjugates."""
    if len(x) != len(y) or x == y:
        return False
    # x = uv, y = vu for non-empty u, v
    # We iterate through all possible split points of x
    for i in range(1, len(x)):
        u = x[:i]
        v = x[i:]
        if y == v + u:
            return True
    return False

def solve():
    """
    Counts the number of finite substitutions (x,y) for words of length <= 3.
    """
    # Step 1: Generate all words of length 0, 1, 2, 3
    words = []
    for i in range(4): # lengths 0, 1, 2, 3
        for j in range(2**i): # all possible words of length i
            # Create binary string and replace 0 with 'a', 1 with 'b'
            words.append(bin(j)[2:].zfill(i).replace('0', 'a').replace('1', 'b'))

    # Step 2 & 3: Count finite substitutions by category
    
    # Case 1: len(y) < len(x)
    # These are always finite.
    finite_len_y_lt_len_x = 0
    for x in words:
        for y in words:
            if len(y) < len(x):
                finite_len_y_lt_len_x += 1
    
    # Case 2: len(y) > len(x)
    # Finite if x is not a subword of y.
    finite_len_y_gt_len_x = 0
    for x in words:
        for y in words:
            if len(y) > len(x):
                # A rule is non-finite if x is a subword of y.
                if not (x in y):
                    finite_len_y_gt_len_x += 1

    # Case 3: len(y) == len(x)
    # Finite if x != y and they are not conjugates.
    # The rule is non-finite if x == y or they are conjugates.
    finite_len_y_eq_len_x = 0
    for x in words:
        for y in words:
            if len(y) == len(x):
                if x != y and not is_conjugate(x, y):
                    finite_len_y_eq_len_x += 1

    total_finite_count = finite_len_y_lt_len_x + finite_len_y_gt_len_x + finite_len_y_eq_len_x
    
    print(f"Total number of possible words for x and y (length <= 3) is {len(words)}.")
    print(f"Total number of substitutions (x, y) is {len(words)} * {len(words)} = {len(words)**2}.")
    print("\nCounting the number of finite substitutions:")
    print(f"1. Finite pairs where len(y) < len(x): {finite_len_y_lt_len_x}")
    print(f"2. Finite pairs where len(y) > len(x): {finite_len_y_gt_len_x}")
    print(f"3. Finite pairs where len(y) == len(x): {finite_len_y_eq_len_x}")
    print("\nThe total number of finite substitutions is the sum of these counts.")
    print(f"Final equation: {finite_len_y_lt_len_x} + {finite_len_y_gt_len_x} + {finite_len_y_eq_len_x} = {total_finite_count}")

solve()
<<<148>>>