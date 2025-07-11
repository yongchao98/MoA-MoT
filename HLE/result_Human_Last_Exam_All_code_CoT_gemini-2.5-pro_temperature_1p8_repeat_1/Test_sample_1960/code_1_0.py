def generate_and_print_solution():
    """
    This function generates and prints the linear logic encoding for the equipartitioning problem.
    """
    
    # 1. Define the base formula U. We choose a simple, literal-free formula.
    U_formula = "(1 ⊸ ⊥)"

    def power_string(base, exp_symbol):
        """Generates a human-readable string for a formula raised to a power."""
        base_repr = base if ' ' not in base.strip() else f"({base})"
        
        # Use a dictionary for special cases for clarity
        cases = {
            "0": "1",
            "1": base,
            "w": f"{base_repr} ⊗ ... ⊗ {base_repr} ({exp_symbol} times)",
            "b": f"{base_repr} ⊗ ... ⊗ {base_repr} ({exp_symbol} times)",
            "m": f"{base_repr} ⊗ ... ⊗ {base_repr} ({exp_symbol} times)"
        }
        
        if exp_symbol in cases:
            return cases[exp_symbol]
        else: # Generic case for any other symbol
            return f"{base_repr} ⊗ ... ⊗ {base_repr} ({exp_symbol} times)"

    print("This script describes the solution for encoding the equipartitioning problem in linear logic.")
    
    # --- Print Base Formula ---
    print("\n" + "="*50)
    print("1. The Base Formula (U)")
    print("-"*50)
    print("We define a base formula U that represents one unit of 'value'.")
    print("This formula must not contain any literals. A suitable choice is:")
    print(f"\n  U = {U_formula}\n")

    # --- Print function f(w) ---
    print("="*50)
    print("2. The Function f(w)")
    print("-"*50)
    print("The function f(w) maps a natural number w to a formula representing w units.")
    print("This is done by tensoring U with itself 'w' times. Each number w from the input set W is converted to a formula f(w).")
    print("For a given number 'w', the formula is:")
    print(f"\n  f(w) = {power_string(U_formula, 'w')}\n")
    print("Note: If w = 0, f(0) is the multiplicative unit '1'. If w = 1, f(1) is just U.")
    
    # --- Print Goal Formula C ---
    print("="*50)
    print("3. The Goal Formula C(W, m, b)")
    print("-"*50)
    print("The goal formula C represents the target state: 'm' partitions, each summing to 'b'.")
    print("A single partition that sums to 'b' is represented by the formula U^b.")
    
    u_b_string = power_string(U_formula, 'b')
    print("\n  Formula for one partition of sum 'b':")
    print(f"  U^b = {u_b_string}\n")
    
    print("The final goal C requires 'm' such partitions simultaneously, so we tensor them together:")
    # The base for this power operation is the string representation of U^b, which is "(U ⊗ ... ⊗ U (b times))"
    c_formula_string = power_string("(U^b)", 'm')
    print("\n  Final goal formula C:")
    print(f"  C(m, b) = (U^b)^m = {c_formula_string}\n")


# Execute the function to print the solution.
generate_and_print_solution()