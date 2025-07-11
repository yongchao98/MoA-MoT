def is_finite(x, y):
    """
    Determines if the substitution x -> y is finite based on a set of rules.
    """
    # Rule a: If x is the empty string, the substitution is not finite.
    if not x:
        return False

    # Rule b: If y is shorter than x, the substitution is finite.
    if len(y) < len(x):
        return True

    # Rule c: If y has the same length as x, it's finite iff x != y.
    if len(y) == len(x):
        return x != y

    # Rule d: If y is longer than x.
    # The substitution is infinite if x is a substring of y.
    if x in y:
        return False

    # Check for infinite loops caused by overlaps.
    # An overlap can occur if a proper part of x matches the beginning or end of y.
    
    # Check if a proper suffix of x is a prefix of y.
    # e.g., x="ab", y="bab". Suffix "b" of x is a prefix of y. Infinite: ab -> bab -> b(ab)
    for i in range(1, len(x)):
        suffix = x[i:]
        if y.startswith(suffix):
            return False

    # Check if a proper prefix of x is a suffix of y.
    # e.g., x="ba", y="abab". Prefix "b" of x is a suffix of y. Infinite: ba -> abab -> a(ba)b
    for i in range(1, len(x)):
        prefix = x[:i]
        if y.endswith(prefix):
            return False

    # If none of the infinite conditions are met, the substitution is finite.
    return True

def solve():
    """
    Counts the number of finite substitutions for words of length <= 3.
    """
    alphabet = ['a', 'b']
    
    # Generate all words of length 0 to 3
    all_words = [""]
    current_gen = [""]
    for _ in range(3):
        next_gen = [word + char for word in current_gen for char in alphabet]
        all_words.extend(next_gen)
        current_gen = next_gen

    # The prompt mentions 255 substitutions, but the described set of words
    # (length <= 3 on {a,b}) generates 15*15 = 225 possible substitutions.
    # We will calculate the answer for these 225 pairs.
    
    total_finite_count = 0
    
    # We group counts by the length of x for a clear breakdown.
    counts_by_len_x = {1: 0, 2: 0, 3: 0}

    for x in all_words:
        if not x:  # The case x="" is never finite and len(x)=0.
            continue
            
        len_x = len(x)
        for y in all_words:
            if is_finite(x, y):
                counts_by_len_x[len_x] += 1
                total_finite_count += 1
    
    print("Analysis of finite substitutions (x -> y):")
    for length, count in counts_by_len_x.items():
        num_x = 2**length
        total_y = len(all_words)
        print(f"For x of length {length} ({num_x} words): there are {count} finite substitutions out of {num_x * total_y}.")

    print(f"\nTotal number of finite substitutions = {counts_by_len_x[1]} + {counts_by_len_x[2]} + {counts_by_len_x[3]}")
    print(f"Total number of finite substitutions: {total_finite_count}")

solve()
<<<150>>>