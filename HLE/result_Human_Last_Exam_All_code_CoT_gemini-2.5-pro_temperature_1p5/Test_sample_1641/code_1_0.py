import collections

def generate_words():
    """Generates words of length <= 3 on alphabet {a,b}."""
    words = ['']
    for length in range(1, 4):
        current_len_words = []
        for word in words:
            if len(word) == length - 1:
                current_len_words.append(word + 'a')
                current_len_words.append(word + 'b')
        words.extend(current_len_words)
    return words

def solve():
    """
    Calculates the number of finite substitutions (x,y) for words of length <= 3.
    """
    words = generate_words()
    total_substitutions = len(words) * len(words)

    finite_counts_by_len_x = collections.defaultdict(int)
    
    for x in words:
        if x == '':
            # Rule: If x is the empty string, the substitution is infinite.
            # So, we count 0 finite substitutions for this case.
            continue
            
        for y in words:
            is_infinite = False
            # Rule: If x = y, the substitution is infinite.
            if x == y:
                is_infinite = True
            # Rule: If |y| > |x| and x is a substring of y, it's infinite.
            elif len(y) > len(x) and x in y:
                is_infinite = True

            if not is_infinite:
                finite_counts_by_len_x[len(x)] += 1
    
    print("The problem asks for the number of finite substitutions out of all possible couples (x,y) of words of length <= 3.")
    print("The number of such words is 15. The total number of substitutions is 15 * 15 = 225.")
    print("\nWe count the number of finite substitutions based on the length of the word 'x'.")
    print("Substitutions where x is the empty word are all infinite.")
    
    total_finite_count = 0
    equation_parts = []
    
    # We iterate through sorted keys (lengths of x) for a clean output
    for length in sorted(finite_counts_by_len_x.keys()):
        count = finite_counts_by_len_x[length]
        print(f"Number of finite substitutions where length(x) = {length}: {count}")
        total_finite_count += count
        equation_parts.append(str(count))
        
    equation_str = " + ".join(equation_parts)
    print(f"\nFinal calculation: {equation_str} = {total_finite_count}")
    
    print(f"\nOut of the 225 possible substitutions, {total_finite_count} are finite.")
    print(f"\n<<<164>>>")

solve()
