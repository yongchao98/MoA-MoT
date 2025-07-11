def represent_formula(base_formula_str: str, n: int) -> str:
    """
    Creates a string representation of a formula tensored with itself n times.
    Example: (F) ⊗ (F) ⊗ ... ⊗ (F)
    """
    if not isinstance(n, int) or n < 0:
        raise ValueError("The number of repetitions (n) must be a non-negative integer.")
    if n == 0:
        # The identity element for the tensor product is 1.
        return "1"
    if n == 1:
        return f"({base_formula_str})"
    
    # For n > 1, join n copies of the base formula with the tensor symbol.
    return " ⊗ ".join([f"({base_formula_str})" for _ in range(n)])

def print_equipartition_encoding(W: list[int], m: int, b: int):
    """
    Prints the linear logic sequent that encodes the specified
    equipartitioning problem EP(W, m, b).
    """
    print(f"Generating the linear logic encoding for the equipartitioning problem:")
    print(f"W = {W}, m = {m}, b = {b}\n")
    
    # Check the necessary condition for a solution to exist.
    if sum(W) != m * b:
        print(f"Error: The sum of elements in W ({sum(W)}) does not equal m * b ({m*b}).")
        print("An equipartition is impossible, and the corresponding sequent is not derivable.")
        return

    # A base formula A constructed from allowed connectives. A = (1 ⊸ ⊥).
    # This formula acts as a unique, non-trivial "counter" resource.
    base_formula = "1 ⊸ ⊥"

    # Define f(w) = A^{\otimes w}
    # These are the initial resources (the antecedent of the sequent).
    f_w_formulas = [represent_formula(base_formula, w) for w in W]
    antecedent = ", ".join(f_w_formulas)

    # Define C = (A^{\otimes b})^{\otimes m}
    # This is the target state (the succedent of the sequent).
    partition_formula = represent_formula(base_formula, b)
    c_formula = represent_formula(f"({partition_formula})", m)
    
    print("---------------------------------------------------------------------")
    print("The provable sequent in linear logic is constructed as follows:")
    print("1. Define a base 'unit' formula: A = (1 ⊸ ⊥)")
    print(f"2. For each number w in W, define a resource formula f(w) = A^{{\otimes w}}")
    print(f"   For w={W[0]}, f({W[0]}) = {f_w_formulas[0]}")
    
    print(f"3. Define the target formula C = (A^{{\otimes b}})^{{\otimes m}}")
    print(f"   For m={m} and b={b}, C = {c_formula}")

    print("\n4. The full sequent is: { f(w) | w ∈ W } ⊢ C")
    print("\nThis expands to:")
    print(f"   {antecedent} \n   ⊢ \n   {c_formula}")
    print("---------------------------------------------------------------------")
    print("\nThis sequent is provable if and only if the set W can be partitioned")
    print(f"into {m} subsets, each summing to {b}.")

# --- Example Usage ---
# A problem instance where a solution exists: W={6, 5, 5, 4}, m=2, b=10
# Partitions are {6, 4} and {5, 5}
W_example = [6, 5, 5, 4]
m_example = 2
b_example = 10
print_equipartition_encoding(W_example, m_example, b_example)