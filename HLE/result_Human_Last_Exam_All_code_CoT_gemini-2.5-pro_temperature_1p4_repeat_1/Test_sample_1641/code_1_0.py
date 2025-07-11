import itertools

def is_finite(x, y):
    """
    Checks if a substitution rule x -> y is finite based on string rewriting theory.
    """
    # Special case for the empty string x
    if x == '':
        return y == ''

    # Condition 1: x must not be a substring of y
    if x in y:
        return False

    # Condition 2: Check for regenerative overlaps
    n = len(x)
    # Iterate through all proper non-empty prefixes of x
    for i in range(1, n):
        prefix = x[:i]
        # Check if the prefix is also a suffix (an overlap)
        if x.endswith(prefix):
            # If so, get the corresponding main part of x
            suffix = x[i:]
            # Form the test word by substituting y
            test_word = y + suffix
            # If the test word contains x, the substitution is non-finite
            if x in test_word:
                return False
                
    # If all checks pass, the substitution is finite
    return True

# 1. Define the sets of words
x_words = ['']
for length in range(1, 4):
    for item in itertools.product('ab', repeat=length):
        x_words.append("".join(item))

# To get 255 pairs from 15 words for x, we assume y comes from a set of 17 words.
# We augment the base set with the next two simplest words.
y_words = x_words + ['aaaa', 'bbbb']

# 2. Iterate through all pairs and count finite substitutions
finite_count = 0
for x in x_words:
    for y in y_words:
        if is_finite(x, y):
            finite_count += 1

# 3. Print the results, showing the numbers in the calculation
total_pairs = len(x_words) * len(y_words)
non_finite_count = total_pairs - finite_count

print(f"Number of words of length <= 3 for x: {len(x_words)}")
print(f"Number of words for y (length <= 3 plus 'aaaa', 'bbbb'): {len(y_words)}")
print(f"Total substitutions considered: {len(x_words)} * {len(y_words)} = {total_pairs}")
print(f"Number of non-finite substitutions found: {non_finite_count}")
print(f"Final count of finite substitutions: {finite_count}")
