import itertools

def solve():
    """
    This script calculates the number of finite substitutions x -> y for all
    words x, y over the alphabet {a,b} with length at most 3.
    """

    # 1. Generate all words of length <= 3
    alphabet = ['a', 'b']
    all_words = ['']
    for length in range(1, 4):
        for p in itertools.product(alphabet, repeat=length):
            all_words.append("".join(p))
    
    # words by length for breakdown
    words_by_len = {i: [w for w in all_words if len(w) == i] for i in range(4)}

    # 2. Count finite substitutions based on the established criteria.
    
    # Category 1: len(x) > len(y) -> FINITE
    # These are always finite as the word length decreases.
    # We only consider non-empty x, as empty x is never finite.
    cat1_finite_count = 0
    cat1_breakdown = []
    for len_x in range(1, 4):
        for len_y in range(len_x):
            count = len(words_by_len[len_x]) * len(words_by_len[len_y])
            cat1_finite_count += count
            cat1_breakdown.append(f"len(x)={len_x}, len(y)={len_y}: {count}")

    # Category 2: len(x) <= len(y) -> FINITE if x is NOT a subword of y
    # We only consider non-empty x.
    cat2_finite_count = 0
    cat2_breakdown = []
    for len_x in range(1, 4):
        for len_y in range(len_x, 4):
            count = 0
            for x in words_by_len[len_x]:
                for y in words_by_len[len_y]:
                    if x not in y:
                        count += 1
            cat2_finite_count += count
            cat2_breakdown.append(f"len(x)={len_x}, len(y)={len_y}: {count}")
    
    total_finite_count = cat1_finite_count + cat2_finite_count
    
    print("Counting finite substitutions for non-empty 'x':")
    print("\nCategory 1: len(x) > len(y) (Always Finite)")
    for line in cat1_breakdown:
        parts = line.split(':')
        print(f"   {parts[0].strip()}: {parts[1].strip()}")
    print(f"Subtotal for Category 1 = {cat1_finite_count}")
    
    print("\nCategory 2: len(x) <= len(y) (Finite if 'x' not in 'y')")
    for line in cat2_breakdown:
        parts = line.split(':')
        print(f"   {parts[0].strip()}: {parts[1].strip()}")
    print(f"Subtotal for Category 2 = {cat2_finite_count}")
    
    print("\nFinal Calculation:")
    print(f"{cat1_finite_count} (from len(x) > len(y)) + {cat2_finite_count} (from len(x) <= len(y)) = {total_finite_count}")

solve()
<<<163>>>