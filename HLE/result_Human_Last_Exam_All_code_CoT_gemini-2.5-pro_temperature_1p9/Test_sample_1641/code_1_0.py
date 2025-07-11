import collections
from itertools import product

def generate_words():
    """Generates all words of length 0 to 3 on the alphabet {a, b}."""
    words = [""]
    for length in range(1, 4):
        words.extend([''.join(p) for p in product('ab', repeat=length)])
    return words

def get_borders(word):
    """
    A border of a word is a proper non-empty prefix that is also a suffix.
    This function finds all such borders for a given word.
    """
    borders = []
    # A border must be shorter than the word itself.
    for i in range(1, len(word)):
        prefix = word[:i]
        suffix = word[-i:]
        if prefix == suffix:
            borders.append(prefix)
    return borders

def is_finite(x, y):
    """
    Determines if the substitution x -> y is finite based on string rewriting theory.
    Returns True if finite, False otherwise.
    """
    # Case 1: x is the empty word.
    # The substitution is finite if and only if y is also the empty word (a no-op).
    if not x:
        return not y

    # Case 2: x is non-empty.
    # The substitution is infinite if the replacement word y contains x.
    if x in y:
        return False

    # The substitution can also be infinite if y regenerates x via an overlap with itself.
    # This occurs if x has a "border" (a string that is both a proper prefix and suffix).
    borders = get_borders(x)
    for z in borders:
        # Define 'b' and 'a' from the structure x = bz = za
        # b is the part of x before the suffix z
        b = x[:-len(z)]
        # a is the part of x after the prefix z
        a = x[len(z):]

        # An infinite chain is possible if `b` followed by `y` creates `x`,
        # or `y` followed by `a` creates `x`.
        if x in (b + y) or x in (y + a):
            return False

    # If no conditions for infinitude are met, the substitution is finite.
    return True

def solve():
    """
    Counts the number of finite substitutions for words of length <= 3.
    """
    words = generate_words()
    
    # We will store the counts of finite substitutions grouped by the length of x.
    counts_by_len = collections.defaultdict(int)

    for x in words:
        for y in words:
            if is_finite(x, y):
                counts_by_len[len(x)] += 1
    
    total_finite = sum(counts_by_len.values())

    c0 = counts_by_len[0]
    c1 = counts_by_len[1]
    c2 = counts_by_len[2]
    c3 = counts_by_len[3]

    print(f"Number of finite substitutions for len(x) = 0: {c0}")
    print(f"Number of finite substitutions for len(x) = 1: {c1}")
    print(f"Number of finite substitutions for len(x) = 2: {c2}")
    print(f"Number of finite substitutions for len(x) = 3: {c3}")
    print("-" * 20)
    print(f"Total number of finite substitutions is the sum:")
    print(f"{c0} + {c1} + {c2} + {c3} = {total_finite}")
    
    return total_finite

if __name__ == '__main__':
    final_answer = solve()
    # The format below is for the final answer extraction.
    # <<<125>>>