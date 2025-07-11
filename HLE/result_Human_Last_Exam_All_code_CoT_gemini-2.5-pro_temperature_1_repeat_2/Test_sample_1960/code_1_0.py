def tensor_power(base, exp, op="⊗"):
    """
    Creates a string for a formula raised to a power with a tensor, e.g., A^3 -> A ⊗ A ⊗ A.
    """
    if exp == 0:
        return "1"
    if exp == 1:
        return base
    return f" {op} ".join([base] * exp)

def solve_equipartition_encoding(W, m, b):
    """
    Generates and prints the linear logic formulas for the equipartitioning problem.
    
    Args:
        W (set): A set of natural numbers.
        m (int): The number of partitions.
        b (int): The target sum for each partition.
    """
    
    # Verify the sum constraint
    if sum(W) != m * b:
        print(f"Error: The sum of elements in W ({sum(W)}) does not equal m * b ({m*b}).")
        print("The equipartitioning problem is ill-defined for these inputs.")
        return

    print("--- Equipartitioning Problem ---")
    print(f"W = {W}")
    print(f"m = {m}")
    print(f"b = {b}\n")
    
    print("--- Linear Logic Encoding ---\n")

    # Define and print the function f(w) for each w in W
    print("The function f(w) maps each number w to a formula representing w resources:")
    for w in sorted(list(W)):
        formula_fw = tensor_power("⊥", w)
        print(f"f({w}) = {formula_fw}")

    print("-" * 20)

    # Define and print the goal formula C
    print("The goal formula C is constructed as follows:")
    
    # 1. Formula for a bin of size b
    bin_b_formula = tensor_power("⊥", b)
    print(f"1. A bin of size {b} is represented by: \n   B = {bin_b_formula}")

    # 2. Formula for a machine that fills one bin
    bin_filler_formula = f"({bin_b_formula} ⊸ 1)"
    print(f"\n2. A 'bin-filler' machine that consumes B and produces 1 is: \n   M = {bin_filler_formula}")

    # 3. Formula for m machines
    m_machines_formula = tensor_power(bin_filler_formula, m)
    if m > 1:
      m_machines_formula = f"({m_machines_formula})"
    print(f"\n3. The complete set of {m} machines is: \n   M_total = {m_machines_formula}")
    
    # 4. Final formula C
    c_formula = f"{m_machines_formula} ⊸ 1"
    print(f"\n4. The final goal formula C is: \n   C = {c_formula}\n")

    print("-" * 20)
    print("The equipartitioning problem EP(W, m, b) is true if and only if the following sequent is derivable in linear logic:")
    
    sequent_lhs = ", ".join([f"f({w})" for w in sorted(list(W))])
    print(f"{{ {sequent_lhs} }} ⊢ C")


# --- Example Usage ---
# You can change these values to see the encoding for different problems.
# Example where the partition exists: W = {1, 2, 3, 4}, m = 2, b = 5
W_example = {1, 2, 3, 4}
m_example = 2
b_example = 5

solve_equipartition_encoding(W_example, m_example, b_example)

<<<
f(w) = ⊥^w (which is ⊥ ⊗ ... ⊗ ⊥, w times; f(0) = 1)
C = ((⊥^b ⊸ 1)^m) ⊸ 1 (where F^k is F ⊗ ... ⊗ F, k times)
>>>