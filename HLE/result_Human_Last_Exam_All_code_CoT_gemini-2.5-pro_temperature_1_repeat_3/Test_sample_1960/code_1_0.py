import argparse

def tensor_power(n, base_formula_str):
    """Generates the n-th tensor power of a formula string."""
    if n < 0:
        raise ValueError("Power must be a non-negative integer.")
    if n == 0:
        return "1"
    if n == 1:
        return base_formula_str
    # Using explicit parentheses for clarity
    return f"({' @ '.join([base_formula_str] * n)})"

def nested_implication(n, hypothesis_str, conclusion_str):
    """Generates a nested implication F --o (F --o ... (F --o G)...)."""
    if n < 0:
        raise ValueError("Nesting level must be a non-negative integer.")
    if n == 0:
        return conclusion_str
    
    res = conclusion_str
    for _ in range(n):
        res = f"({hypothesis_str} --o {res})"
    return res

def generate_formulas(W, m, b):
    """
    Prints the string representation of the linear logic formulas
    for the equipartitioning problem.
    """
    # Define a base formula 'A' that contains no literals.
    # A = (bot --o 1) is a suitable choice as it's not equivalent to 1 or bot.
    A = "(bot --o 1)"

    print("--- Encoding Definition ---")
    print(f"Base formula A = {A}")
    print("Function f(w) = A^w (w-th tensor power of A)")
    print("Formula C = P --o (P --o ... (P --o 1)...) (m times)")
    print("where P = A^b --o 1 (b-th tensor power of A)")
    print("-" * 25)
    print("\n--- Generated Formulas ---")
    
    # Print f(w) for each w in W
    print("Formulas for the context {f(w) | w in W}:")
    if not W:
        print("W is empty.")
    else:
        for w in W:
            f_w = tensor_power(w, A)
            print(f"f({w}) = {f_w}")
    
    print("\nFormula for the conclusion C:")
    # Define the formula for a partition of size b
    B = tensor_power(b, A)
    P = f"({B} --o 1)"
    
    # Define the final formula C
    C = nested_implication(m, P, "1")
    
    print(f"For m = {m}, b = {b}:")
    print(f"C = {C}")


if __name__ == '__main__':
    # We don't need to actually solve the problem, just show the formulas.
    # So we use some example values for W, m, b.
    # Example: W = {2, 4}, m = 2, b = 3. sum(W) = 6, m*b = 6.
    example_W = [2, 4]
    example_m = 2
    example_b = 3
    
    # Check the problem constraint
    if sum(example_W) != example_m * example_b:
        print(f"Warning: The condition sum(W) = m*b is not met for the example values.")
        print(f"sum({example_W}) = {sum(example_W)}")
        print(f"{example_m} * {example_b} = {example_m * example_b}")
    else:
        generate_formulas(example_W, example_m, example_b)
