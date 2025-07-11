import collections

def generate_formula_string(base_formula, count):
    """
    Generates a string for (A ⊗ A ⊗ ... ⊗ A) with proper parenthesizing.
    """
    if count == 0:
        return "1"
    if count == 1:
        return base_formula
    
    # Using the Unicode character for the tensor product (⊗)
    return " ⊗ ".join([base_formula] * count)

def generate_ep_sequent(W, m, b):
    """
    Prints the linear logic sequent corresponding to the equipartitioning problem EP(W, m, b).
    This function outputs each number from the inputs W, m, and b in the final sequent representation.
    """
    # Verify the sum constraint for the problem's validity
    if sum(W) != m * b:
        print(f"Error: The sum of elements in W ({sum(W)}) does not equal m*b ({m*b}).")
        print("The equipartitioning problem is not well-defined for these inputs.")
        return

    # 1. Define the base formula 'A'. Using Unicode for logic symbols:
    # ⊸ (lollipop), ⊥ (bottom), ⊗ (tensor), ⊢ (turnstile)
    A = "(1 ⊸ ⊥)"

    # 2. Generate the set of antecedent formulas, Γ = {f(w) | w ∈ W}
    # Using Counter to handle duplicate numbers in W correctly
    w_counts = collections.Counter(W)
    antecedent_formulas_str = []
    for w, count in sorted(w_counts.items()):
        # Generate f(w) = A^w
        f_w = generate_formula_string(A, w)
        # Add the formula 'count' times to our list of resources
        antecedent_formulas_str.extend([f"f({w}) = {f_w}"] * count)

    # 3. Generate the succedent formula, C = (A^b)^m
    # Formula for one bucket summing to 'b'
    B = generate_formula_string(A, b)
    
    # Formula for 'm' such buckets
    # Parenthesize B if it's a tensor product itself, for clarity
    B_paren = f"({B})" if b > 1 else B
    C = generate_formula_string(B_paren, m)

    # 4. Print the final, formatted sequent
    print("--------------------------------------------------")
    print(f"Encoding for EP(W={W}, m={m}, b={b})")
    print("--------------------------------------------------")
    print("Base Formula (A):")
    print(f"  {A}\n")
    
    print("Function f(w):")
    print(f"  f(w) = {A} ⊗ ... ⊗ {A} (w times)\n")
    
    print("Antecedent (Γ - Resources):")
    print("  {")
    for formula_str in antecedent_formulas_str:
        print(f"    {formula_str}")
    print("  }\n")
    
    print("Succedent (C - Goal):")
    print(f"  {C}\n")

    print("Full Sequent (Γ ⊢ C):")
    # Join antecedent formulas with ', ' for the final sequent string
    antecedent_str = ",\n ".join([f"({f.split('= ')[1]})" for f in antecedent_formulas_str])
    print(f"{antecedent_str}")
    print(f"⊢")
    print(f"{C}")
    print("--------------------------------------------------")

# Example Usage:
# Define a solvable equipartitioning problem instance
W = {6, 5, 5, 4} # sum(W) = 20
m = 2             # number of partitions
b = 10            # target sum for each partition
# Note: m * b = 2 * 10 = 20, which matches sum(W).
# A valid partition is W1={6, 4} and W2={5, 5}.

generate_ep_sequent(W, m, b)
