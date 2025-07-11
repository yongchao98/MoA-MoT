def f_to_string(w):
    """
    Generates the string representation for the formula f(w).
    f(w) is defined as bot tensored w times.
    """
    if not isinstance(w, int) or w < 0:
        raise ValueError("w must be a non-negative integer.")
    if w == 0:
        return "1"
    return " * ".join(["bot"] * w)

def C_to_string(m, b):
    """
    Generates the string representation for the formula C(W, m, b).
    C is B tensored m times, where B checks for a sum of b.
    """
    if not isinstance(m, int) or m < 0:
        raise ValueError("m must be a non-negative integer.")
    if not isinstance(b, int) or b < 0:
        raise ValueError("b must be a non-negative integer.")

    # Base case: C is 1 if m is 0.
    if m == 0:
        return "1"

    # 1. Construct the formula for a single bin, B.
    # B = ((bot^{\otimes b} -> bot) -> bot)
    
    # Inner part: bot^{\otimes b}
    if b == 0:
        bot_b_str = "1"
    else:
        # This part explicitly shows the number 'b' in the formula
        bot_b_str = " * ".join(["bot"] * b)

    B_str = f"(( {bot_b_str} -> bot ) -> bot)"

    # 2. Construct the final formula C by tensoring B m times.
    # This part explicitly shows the number 'm' in the formula
    C_str = " * ".join([B_str] * m)
    
    return C_str

def print_encoding_definitions(m, b):
    """
    Prints the definitions for f(w) and C(m,b) in a human-readable format.
    We use 'bot' for the constant ⊥, '*' for ⊗, and '->' for ⊸.
    """
    print("--- Linear Logic Encoding for Equipartitioning ---")
    print("\n1. The function f(w):")
    print("f(w) maps a natural number w to w copies of the 'bot' constant, combined by the tensor connective '*'.")
    print("f(w) = bot * bot * ... * bot (w times)")
    print("f(0) = 1")
    print("\nExample for w=3:")
    print(f"f(3) = {f_to_string(3)}")
    
    print("\n2. The formula C(m, b):")
    print("C is the goal formula, representing m 'bins', each needing a sum of exactly b.")
    print(f"\nFor the specific instance with m = {m} and b = {b}, the formula C is:")
    print(f"C = {C_to_string(m, b)}")
    print("\nEach parenthesized block `((...) -> bot)` represents one bin.")
    print(f"The formula contains {m} such blocks, corresponding to m partitions.")
    print(f"Inside each block, there are {b} instances of 'bot', corresponding to the target sum b.")


# Example usage:
# You can change these values to see the formulas for a different problem instance.
m_example = 2
b_example = 4

print_encoding_definitions(m_example, b_example)
