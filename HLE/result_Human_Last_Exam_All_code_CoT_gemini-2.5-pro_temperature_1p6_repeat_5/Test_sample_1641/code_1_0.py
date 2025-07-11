import itertools

def count_finite_substitutions():
    """
    Counts the number of finite substitutions x -> y for words x, y
    of length <= 3 on the alphabet {a, b}.
    """
    
    # 1. Generate all words of length <= 3
    words = [""]
    for length in range(1, 4):
        words.extend([''.join(p) for p in itertools.product('ab', repeat=length)])

    # 2. Initialize counters for each case
    finite_gt = 0  # for len(x) > len(y)
    finite_lt = 0  # for len(x) < len(y)
    finite_eq = 0  # for len(x) == len(y)

    # 3. Iterate through all (x, y) pairs and apply finiteness rules
    for x in words:
        for y in words:
            is_finite = False
            
            # Rule for x = "" (always non-finite)
            if not x:
                is_finite = False
            # Rule for len(x) > len(y)
            elif len(x) > len(y):
                is_finite = True
            # Rule for len(x) < len(y)
            elif len(x) < len(y):
                if x not in y:
                    is_finite = True
            # Rule for len(x) == len(y)
            else: # len(x) == len(y) and x != ""
                if y < x: # Lexicographical comparison
                    is_finite = True

            if is_finite:
                if len(x) > len(y):
                    finite_gt += 1
                elif len(x) < len(y):
                    finite_lt += 1
                else: # len(x) == len(y)
                    finite_eq += 1
    
    # 4. Output the results in the form of an equation
    total_finite = finite_gt + finite_lt + finite_eq
    print(f"Number of finite substitutions where:")
    print(f"len(x) > len(y): {finite_gt}")
    print(f"len(x) < len(y): {finite_lt}")
    print(f"len(x) = len(y): {finite_eq}")
    print("---")
    print(f"Total: {finite_gt} + {finite_lt} + {finite_eq} = {total_finite}")

count_finite_substitutions()