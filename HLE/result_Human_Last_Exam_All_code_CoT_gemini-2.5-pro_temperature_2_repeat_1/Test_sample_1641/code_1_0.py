import collections

def main():
    """
    Calculates the number of finite substitutions x -> y for words x, y
    on the alphabet {a,b} with length <= 3.
    """

    # 1. Generate all words of length <= 3
    max_len = 3
    alphabet = ['a', 'b']
    words = [""]
    q = collections.deque([""])
    # Generate words breadth-first
    while q:
        curr = q.popleft()
        if len(curr) < max_len:
            for char in alphabet:
                new_word = curr + char
                words.append(new_word)
                q.append(new_word)

    # 2. Define helper function for finiteness conditions
    def is_conjugate(s1, s2):
        """Checks if s2 is a cyclic permutation of s1."""
        if not s1 or len(s1) != len(s2):
            return False
        return s2 in s1 + s1

    # 3. Iterate through all pairs and count finite substitutions
    count_len_decrease = 0
    count_len_increase_finite = 0
    count_len_equal_finite = 0

    for x in words:
        for y in words:
            len_x = len(x)
            len_y = len(y)

            # Case: |y| < |x| => Finite
            if len_y < len_x:
                count_len_decrease += 1
                continue

            # Case: x is empty word => Infinite
            if not x: # x is ""
                continue
            
            # Now, |y| >= |x| and x is non-empty

            # Case: |y| > |x| => Finite if x is not a substring of y
            if len_y > len_x:
                if x not in y:
                    count_len_increase_finite += 1
            
            # Case: |y| = |x| => Finite if y is not a conjugate of x
            elif len_y == len_x:
                if not is_conjugate(x, y):
                    count_len_equal_finite += 1
    
    total_finite = count_len_decrease + count_len_increase_finite + count_len_equal_finite
    
    print(f"{count_len_decrease} + {count_len_increase_finite} + {count_len_equal_finite} = {total_finite}")

if __name__ == "__main__":
    main()