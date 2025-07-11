def generate_counterexample_strings(max_i, max_k):
    """
    Generates strings for the language L' = {a^i b^(i+k) a^k}
    which is a subset of the conjugates of K = {a^n b^n}.
    This language L' is not context-free.

    Args:
        max_i (int): The maximum value for the exponent i.
        max_k (int): The maximum value for the exponent k.
    """
    print(f"Generating strings for the non-context-free language L' = {{a^i b^(i+k) a^k}}")
    print("This demonstrates the counterexample for question (c).")
    print("-" * 60)
    
    examples = []
    for i in range(1, max_i + 1):
        for k in range(1, max_k + 1):
            num_b = i + k
            word = 'a' * i + 'b' * num_b + 'a' * k
            examples.append((word, i, k, num_b))

    for word, i, k, num_b in examples:
        # The prompt requires outputting the numbers in the final equation.
        # The equation for the exponents is: num_b = i + k
        print(f"Generated word: \"{word}\"")
        print(f"Equation: {num_b} = {i} + {k}")
        print("-" * 20)

if __name__ == '__main__':
    # Generate a few examples
    generate_counterexample_strings(max_i=3, max_k=2)
