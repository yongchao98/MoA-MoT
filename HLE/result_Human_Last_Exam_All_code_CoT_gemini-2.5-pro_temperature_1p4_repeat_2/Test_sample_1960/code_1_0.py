def tensor_product(base_formula: str, n: int) -> str:
    """Creates a tensor product of a formula n times."""
    if n < 0:
        raise ValueError("Cannot have a negative number of tensor products.")
    if n == 0:
        return "1"
    if n == 1:
        return base_formula
    
    # For n > 1, wrap each element in parentheses for clarity
    elements = [f"({base_formula})" for _ in range(n)]
    return " @ ".join(elements)

def main():
    """
    Generates and prints the linear logic formulas for the equipartition problem.
    """
    # Use --o for linear implication (lollipop) and @ for tensor product.
    A_formula = "bot"
    B_formula = "(bot --o bot)"

    print("--- General Formula Definitions ---")
    print("\nLet A and B be two orthogonal formulas defined as:")
    print(f"  A = {A_formula}")
    print(f"  B = {B_formula}")
    print("\nWhere '--o' represents linear implication (lollipop, \multimap)")
    print("and '@' represents the tensor product (\otimes).")

    print("\nThe function f(w) maps a natural number w to the formula:")
    print("  f(w) = (B^w) --o (A^w)")
    print("  where X^n is the tensor product of X with itself n times (e.g., B @ B @ ... @ B).")

    print("\nThe goal formula C(W, m, b) is defined as:")
    print("  C = G^m")
    print("  where G is the sub-goal for a single partition, G = (B^b) --o (A^b).")
    
    # --- Example Instantiation ---
    W = {2, 4, 3, 3}
    m = 2
    b = 6
    
    # Verify the sum condition
    if sum(W) != m * b:
        print("\n--- Example: Invalid Problem ---")
        print(f"The given W={W}, m={m}, b={b} is not a valid equipartition problem,")
        print(f"because sum(W) = {sum(W)} != m*b = {m*b}.")
        return

    print("\n\n--- Example: W = {2, 4, 3, 3}, m = 2, b = 6 ---")
    
    print("\nThe required f(w) formulas for each w in W are:")
    for w in sorted(list(W)):
        B_w = tensor_product(B_formula, w)
        A_w = tensor_product(A_formula, w)
        f_w_formula = f"({B_w}) --o ({A_w})"
        print(f"  f({w}) = {f_w_formula}")
        
    print(f"\nThe goal formula C for m={m}, b={b} is:")
    
    B_b = tensor_product(B_formula, b)
    A_b = tensor_product(A_formula, b)
    G_formula = f"(({B_b}) --o ({A_b}))"
    
    C_formula = tensor_product(G_formula, m)
    print(f"  C = {C_formula}")

if __name__ == "__main__":
    main()