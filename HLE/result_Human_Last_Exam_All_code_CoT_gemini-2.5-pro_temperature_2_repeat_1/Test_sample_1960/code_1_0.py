def generate_formula_string(base, power_variable):
    """Generates a human-readable string for a tensor power."""
    if power_variable == '0':
        return "1 (the unit for ⊗)"
    if power_variable == '1':
        return base
    return f"{base} ⊗ ... ⊗ {base} ({power_variable} times)"

def print_solution():
    """Prints the linear logic encoding for the equipartitioning problem."""
    
    # Define the formula for f(w)
    f_w_formula = generate_formula_string("⊥", "w")
    
    # Define the formula for a single bin B
    B_formula = generate_formula_string("⊥", "b")
    
    # Define the formula for C
    C_formula = generate_formula_string("B", "m")

    print("The solution involves defining a function f(w) and a formula C as follows:")
    print("-" * 75)
    
    # Print f(w)
    print("Function f(w):")
    print(f"  f(w) = {f_w_formula}")
    print("\n  This function maps a natural number w to the w-th tensor power of the")
    print("  multiplicative constant ⊥ (bottom). Each number in the set W becomes a")
    print("  resource corresponding to its value.")
    
    print("-" * 75)
    
    # Print C
    print("Formula C = C(W, m, b):")
    print(f"  C = {C_formula}")
    print(f"      where the formula B represents a single bin of size b:")
    print(f"      B = {B_formula}")
    print("\n  The formula C represents the target configuration: m bins, each requiring")
    print("  resources summing to b.")
    
    print("-" * 75)
    
    # Print justification
    print("Justification:")
    print("The sequent {f(w) | w ∈ W} ⊢ C is provable if and only if EP(W, m, b) is true.\n")
    print("1. Provability => EP is true:")
    print("   - To prove C = B ⊗ ... ⊗ B (m times), the proof rules for ⊗ require")
    print("     partitioning the resources {f(w)} into m multisets: Γ₁, ..., Γₘ.")
    print("     This induces a partition of the set W into W₁, ..., Wₘ.")
    print("   - For each partition k, the sequent Γₖ ⊢ B must be provable.")
    print("   - The resources in Γₖ are {⊥^⊗w | w ∈ Wₖ}, which combine to ⊥^⊗(Σw).")
    print("   - The goal B is ⊥^⊗b. The sequent is ⊥^⊗(Σw) ⊢ ⊥^⊗b.")
    print("   - In linear logic, this is only provable if the resources match exactly, so Σw = b.")
    print("   - This must hold for all m partitions, which is the definition of EP(W, m, b).")
    
    print("\n2. EP is true => Provability:")
    print("   - If EP(W, m, b) is true, a partition W₁, ..., Wₘ exists where each ΣWₖ = b.")
    print("   - We can construct the proof by partitioning the resources {f(w)} according to")
    print("     the known partition of W. For each subset Γₖ, the sequent ⊥^⊗b ⊢ ⊥^⊗b is")
    print("     an axiom, so it is provable. Combining the m proofs gives a proof of C.")

if __name__ == '__main__':
    print_solution()