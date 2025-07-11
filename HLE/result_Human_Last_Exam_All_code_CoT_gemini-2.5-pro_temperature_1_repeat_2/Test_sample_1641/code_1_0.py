def solve():
    """
    Counts the number of finite substitutions (x,y) for words on {a,b} of length <= 3.
    """

    def generate_words(alphabet, max_len):
        """Generates all words up to a specified maximum length."""
        words = {""}
        for i in range(max_len + 1):
            words.update({word + char for word in words if len(word) == i -1 for char in alphabet})
        
        # A more straightforward way to generate the words:
        all_words = [""]
        queue = [""]
        while queue:
            word = queue.pop(0)
            if len(word) < max_len:
                for char in alphabet:
                    new_word = word + char
                    all_words.append(new_word)
                    queue.append(new_word)
        return all_words

    def is_finite(x, y):
        """
        Determines if the substitution x -> y is finite based on analytical conditions.
        """
        # If x is the empty word, the substitution can be applied infinitely.
        if not x:
            return False

        lx, ly = len(x), len(y)

        # If the substitution makes the word shorter, it must terminate. Finite.
        if lx > ly:
            return True

        # If lengths are equal, it terminates iff y is different from x.
        if lx == ly:
            return x != y

        # Case: lx < ly. The word can grow.
        # It's infinite if a new 'x' can be created (regeneration).

        # 1. Regeneration by `y` alone: `y` contains `x`.
        if x in y:
            return False

        # 2. Regeneration by context: prefix + y.
        # Check if `p+y` contains `x` for any proper prefix `p` of `x`.
        for i in range(1, lx):
            prefix = x[:i]
            if x in (prefix + y):
                return False

        # 3. Regeneration by context: y + suffix.
        # Check if `y+s` contains `x` for any proper suffix `s` of `x`.
        for i in range(1, lx):
            suffix = x[i:]
            if x in (y + suffix):
                return False

        # If no infinite conditions are met, it's finite.
        return True

    alphabet = ['a', 'b']
    max_len = 3
    words = generate_words(alphabet, max_len)
    
    total_pairs = len(words) * len(words)
    finite_count = 0

    for x in words:
        for y in words:
            if is_finite(x, y):
                finite_count += 1
    
    print(f"Alphabet: {alphabet}")
    print(f"Maximum word length: {max_len}")
    print(f"Total number of words of length <= {max_len}: {len(words)}")
    print(f"Total number of substitutions (x,y) to check: {total_pairs}")
    print(f"Number of finite substitutions found: {finite_count}")

solve()
<<<187>>>