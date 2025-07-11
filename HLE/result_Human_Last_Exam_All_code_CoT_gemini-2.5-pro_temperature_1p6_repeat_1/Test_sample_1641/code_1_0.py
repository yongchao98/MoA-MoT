import itertools

def solve_and_explain():
    """
    Calculates the number of finite substitutions for words of length <= 3
    on the alphabet {a, b} and explains the reasoning.
    """
    alphabet = ['a', 'b']
    max_len = 3

    # Step 1: Generate all words of length <= 3
    words = {""}
    for length in range(1, max_len + 1):
        for p in itertools.product(alphabet, repeat=length):
            words.add("".join(p))
    
    # Sort for deterministic processing
    sorted_words = sorted(list(words), key=lambda w: (len(w), w))
    
    # Step 2: Determine total pairs
    total_words = len(sorted_words)
    total_pairs = total_words * total_words

    print(f"First, we identify all words on {{a, b}} with length from 0 to 3.")
    print(f"Number of words: 1 (len 0) + 2 (len 1) + 4 (len 2) + 8 (len 3) = {total_words}")
    print(f"Total possible pairs (x, y) is {total_words} * {total_words} = {total_pairs}.")
    print("(Note: The prompt's 255 appears to be a typo for 225).\n")

    print("A substitution x -> y is finite if x is not a proper subword of y.")
    print("We will count the number of *infinite* pairs and subtract from the total.\n")

    infinite_counts = []
    print("Counting infinite pairs (where x is a proper subword of y):")
    for x in sorted_words:
        # For a given x, count how many y in our set contain x as a proper subword
        current_infinite_count = 0
        for y in sorted_words:
            if x != y and x in y:
                current_infinite_count += 1
        
        if current_infinite_count > 0:
            x_repr = f"'{x}'" if x else "'' (the empty word)"
            print(f"- For x = {x_repr:<19}, there are {current_infinite_count} words y that contain it as a proper subword.")
            infinite_counts.append(current_infinite_count)

    total_infinite = sum(infinite_counts)
    
    # Step 5: Final calculation output
    equation_str = " + ".join(map(str, infinite_counts))
    print(f"\nTotal number of infinite substitutions = {equation_str} = {total_infinite}")

    final_answer = total_pairs - total_infinite
    print(f"The number of finite substitutions is the total pairs minus the infinite pairs.")
    print(f"Final calculation: {total_pairs} - {total_infinite} = {final_answer}")
    
    return final_answer

if __name__ == "__main__":
    answer = solve_and_explain()
    print(f"\n<<< {answer} >>>")