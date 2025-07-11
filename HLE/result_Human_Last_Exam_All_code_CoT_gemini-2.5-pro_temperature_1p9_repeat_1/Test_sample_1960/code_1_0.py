def generate_tensor_power(base: str, n: int, op: str) -> str:
    """
    Generates a string 'base op base op ... op base' n times.
    Example: generate_tensor_power("bot", 3, "@") -> "bot @ bot @ bot"
    """
    if n < 0:
        raise ValueError("Power must be non-negative.")
    if n == 0:
        # The identity for tensor is 1
        return "1"
    if n == 1:
        return base
    return f" {op} ".join([base] * n)

def f(w: int) -> str:
    """
    Generates the linear logic formula f(w) for a natural number w.
    f(w) = bot ⊗ bot ⊗ ... ⊗ bot (w times)
    """
    # Using '->' for multimap (linear implication) and '@' for otimes (tensor)
    return generate_tensor_power("bot", w, "@")

def C(m: int, b: int) -> str:
    """
    Generates the linear logic formula C for m partitions of sum b.
    C = (P -> (P -> ... -> 1)) where P = (bot^b -> 1)
    """
    if m < 0 or b < 0:
        raise ValueError("m and b must be non-negative.")

    bot_power_b = generate_tensor_power("bot", b, "@")
    
    # Add parentheses for clarity if b > 1
    if b > 1:
        bot_power_b = f"({bot_power_b})"
    
    partition_formula = f"({bot_power_b} -> 1)"

    # Build the nested implication C = P -> (P -> ... (P -> 1)...)
    result = "1"
    if m == 0:
      # If m=0, the problem is trivially true (sum(W) must be 0) and C is just 1.
      return "1"

    for _ in range(m):
        result = f"({partition_formula} -> {result})"
        
    return result

if __name__ == '__main__':
    # Example from the equipartitioning problem.
    # W can be partitioned into {3} and {1, 2}, both summing to b=3.
    W = {1, 2, 3}
    m = 2
    b = 3

    print("Solving for the encoding of the equipartitioning problem instance:")
    print(f"  W = {W}")
    print(f"  m = {m}")
    print(f"  b = {b}")
    print("-" * 20)
    print(f"Checking condition: sum(W) = {sum(W)}, m * b = {m * b}")
    assert sum(W) == m * b, "Sum of elements in W must equal m * b"
    print("-" * 20)

    print("The formulas f(w) for each w in W are:")
    for w in sorted(list(W)):
        formula_f = f(w)
        print(f"f({w}) = {formula_f}")
    
    print("\nThe goal formula C(m, b) is:")
    formula_C = C(m, b)
    print(f"C = {C}")