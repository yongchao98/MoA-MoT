import itertools

def generate_words(alphabet, max_len):
    """Generates all words on the alphabet up to a maximum length."""
    words = []
    for length in range(max_len + 1):
        if length == 0:
            words.append('')
            continue
        # Create all combinations of the given length
        for p in itertools.product(alphabet, repeat=length):
            words.append("".join(p))
    return words

def is_cyclic_shift(x, y):
    """Checks if y is a non-trivial cyclic shift of x."""
    # Must be same length, non-identical, and y must be a substring of x concatenated with itself.
    return len(x) == len(y) and x != y and y in x + x

def is_finite(x, y):
    """
    Determines if the substitution x -> y is finite based on termination rules
    for a one-rule string rewriting system.
    """
    # Case 1: x is the empty word. Finite only if y is also empty.
    if not x:
        return not y

    len_x = len(x)
    len_y = len(y)

    # Case 2: The word gets shorter. Always finite.
    if len_y < len_x:
        return True

    # Case 3: The word gets longer. Finite iff x is not a substring of y.
    if len_y > len_x:
        return x not in y

    # Case 4: The word length is constant.
    if len_y == len_x:
        # If x -> x, it's finite (it does nothing).
        if x == y:
            return True
        # Otherwise, finite iff y is not a non-trivial cyclic shift of x.
        return not is_cyclic_shift(x, y)
    
    # This part of the code should not be reachable.
    return False

def solve():
    """
    Calculates the number of finite substitutions for words of length <= 3
    on the alphabet {a, b}.
    """
    alphabet = ['a', 'b']
    max_length = 3
    words = generate_words(alphabet, max_length)

    # We will store the counts for each length of x in a dictionary.
    counts_by_len_x = {}
    
    for len_x_val in range(max_length + 1):
        # Filter words for the current length of x
        words_x = [w for w in words if len(w) == len_x_val]
        count = 0
        for x in words_x:
            for y in words:
                if is_finite(x, y):
                    count += 1
        counts_by_len_x[len_x_val] = count

    # Get the counts for the final equation
    c0 = counts_by_len_x[0]  # x = ''
    c1 = counts_by_len_x[1]  # len(x) = 1
    c2 = counts_by_len_x[2]  # len(x) = 2
    c3 = counts_by_len_x[3]  # len(x) = 3
    total_finite = c0 + c1 + c2 + c3

    # Print the breakdown and the final result
    print(f"{c0} + {c1} + {c2} + {c3} = {total_finite}")

solve()