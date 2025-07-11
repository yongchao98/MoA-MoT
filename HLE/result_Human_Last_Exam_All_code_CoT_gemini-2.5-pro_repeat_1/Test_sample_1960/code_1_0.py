def tensor_n(formula, n, one_unit='1'):
    """Creates a tensor product of a formula n times."""
    if n < 0:
        raise ValueError("n must be non-negative")
    if n == 0:
        return one_unit
    if n == 1:
        return formula
    return "({})".format(" @ ".join([formula] * n))

def generate_formulas(W, m, b):
    """
    Generates and prints the linear logic formulas for the equipartitioning problem.
    Uses '@' for tensor (⊗) and '-o' for linear implication (⊸).
    """
    # 1. Define the literal-free base formulas
    # O = 1 -o bot
    # A = O (counting unit)
    # B = O -o O (partition catalyst)
    # K = (O -o O) -o O (completed partition token)
    base_O = "(1 -o bot)"
    formula_A = base_O
    formula_B = f"({base_O} -o {base_O})"
    formula_K = f"({formula_B} -o {base_O})"

    print("--- Base Formulas ---")
    print(f"A (counting unit) := {formula_A}")
    print(f"B (partition catalyst) := {formula_B}")
    print(f"K (completed partition) := {formula_K}")
    print("-" * 22)
    print("")

    # 2. Define and print the function f(w) for each w in W
    print("--- Formulas f(w) for w in W ---")
    for w in W:
        A_w = tensor_n(formula_A, w)
        f_w = f"({formula_B} -o ({formula_B} @ {A_w}))"
        print(f"f({w}) := {f_w}")
    print("-" * 34)
    print("")

    # 3. Define and print the goal formula C(m, b)
    print("--- Goal Formula C(m, b) ---")
    A_b = tensor_n(formula_A, b)
    # Machine M = (B @ A^b) -o K
    formula_M = f"(({formula_B} @ {A_b}) -o {formula_K})"
    
    # B^m
    B_m = tensor_n(formula_B, m)
    # M^m
    M_m = tensor_n(formula_M, m)
    # K^m
    K_m = tensor_n(formula_K, m)

    # C = (B^m @ M^m) -o K^m
    formula_C = f"(({B_m} @ {M_m}) -o {K_m})"
    print(f"C({m}, {b}) := {formula_C}")
    print("-" * 28)

if __name__ == '__main__':
    # Example instance of the equipartitioning problem
    # W = {1, 2, 3, 4}, m = 2, b = 5
    # Valid partitions: {1, 4} and {2, 3}
    W_example = [1, 2, 3, 4]
    m_example = 2
    b_example = 5
    
    # Check if the sum constraint holds
    if sum(W_example) != m_example * b_example:
        print("Error: The sum of elements in W must be equal to m * b.")
    else:
        print(f"Encoding EP(W={W_example}, m={m_example}, b={b_example}) into Linear Logic:")
        print("Connectives: '@' is ⊗ (tensor), '-o' is ⊸ (linear implication)\n")
        generate_formulas(W_example, m_example, b_example)
