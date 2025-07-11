import math

# --- User-defined example ---
# EP(W, m, b): Can W be partitioned into m subsets each summing to b?
# Example: W = {1, 2, 3}, m = 2, b = 3.
# sum(W) = 6, m*b = 6.
# This is TRUE, as W can be partitioned into {1, 2} and {3}.
W = {1, 2, 3}
m = 2
b = 3

# --- Formula Generation Logic ---

# Use more readable symbols for the connectives in the output string
TENSOR = "⊗"
IMPLIES = "⊸"
BOT = "⊥"
ONE = "1"

# Memoization caches to avoid recomputing complex formulas
p_formula_cache = {}

def get_binary_representation(n, num_bits):
    """Returns the binary string of n with a fixed number of bits."""
    return format(n, '0' + str(num_bits) + 'b')

def build_p(i, b):
    """
    Builds the formula P_i, which represents the integer i.
    The formula is based on the binary representation of i.
    """
    if (i, b) in p_formula_cache:
        return p_formula_cache[(i, b)]

    if b < 0:
        raise ValueError("b must be non-negative")
    if i < 0 or i > b:
        raise ValueError(f"i must be between 0 and {b}")

    # The two "base bits" for our binary encoding. They are orthogonal.
    base_formula_0 = BOT
    base_formula_1 = f"({BOT} {IMPLIES} {BOT})"

    # Determine the number of bits needed to represent numbers up to b
    if b == 0:
        num_bits = 1
    else:
        num_bits = math.ceil(math.log2(b + 1))
        
    binary_str = get_binary_representation(i, num_bits)
    
    # We build the formula from the inside out, based on the bits
    # P_i = H_{c_0}(H_{c_1}(...(H_{c_{d-1}}(1))...))
    # where c_0 is the least significant bit.
    
    current_formula = ONE
    for bit in reversed(binary_str):
        if bit == '0':
            base = base_formula_0
        else:
            base = base_formula_1
        current_formula = f"({base} {IMPLIES} {current_formula})"
        
    p_formula_cache[(i, b)] = current_formula
    return current_formula

def build_f(w, b):
    """
    Builds the formula f(w).
    f(w) is the tensor product of (P_i --o P_{i+w}) for all possible i.
    """
    if w <= 0:
        raise ValueError("Elements of W must be positive natural numbers.")
        
    # Collect all the implication formulas
    implications = []
    for i in range(b - w + 1):
        p_i = build_p(i, b)
        p_i_plus_w = build_p(i + w, b)
        implications.append(f"({p_i} {IMPLIES} {p_i_plus_w})")
    
    # Join them with the tensor connective
    if not implications:
        return ONE # If w > b, f(w) is just the identity
    return f" {TENSOR} ".join(f"({imp})" for imp in implications)

def build_c(m, b):
    """
    Builds the goal formula C.
    C is the tensor product of m copies of the goal for a single partition.
    """
    if m <= 0:
        raise ValueError("m must be a positive natural number.")
        
    # The goal for one partition is to transform state 0 to state b
    p_0 = build_p(0, b)
    p_b = build_p(b, b)
    single_goal = f"({p_0} {IMPLIES} {p_b})"
    
    # We need to achieve this goal m times
    all_goals = [f"({single_goal})" for _ in range(m)]
    
    return f" {TENSOR} ".join(all_goals)

def main():
    """
    Main function to generate and print the sequent.
    """
    print(f"Encoding the equipartitioning problem EP(W={W}, m={m}, b={b}) in linear logic.")
    print("-" * 20)
    
    # Verify the condition sum(W) = m * b
    if sum(W) != m * b:
        print(f"Error: The condition sum(W) = m*b is not met.")
        print(f"sum({W}) = {sum(W)}, but m*b = {m*b}")
        return

    # --- Left-hand side of the sequent ---
    print("Left-hand side resources {f(w) | w in W}:\n")
    lhs_formulas = []
    for w in sorted(list(W)):
        f_w = build_f(w, b)
        print(f"f({w}) = {f_w}\n")
        lhs_formulas.append(f"f({w})")

    # --- Right-hand side of the sequent ---
    print("-" * 20)
    print("Right-hand side goal C:\n")
    c_formula = build_c(m, b)
    print(f"C = {c_formula}\n")
    
    # --- The full sequent ---
    print("-" * 20)
    print("The complete sequent is:\n")
    lhs_str = ", ".join(lhs_formulas)
    print(f"{lhs_str}  ⊢  C")
    print("\nThis sequent is derivable if and only if the equipartitioning problem is true.")
    print("-" * 20)
    print("Here are the definitions of the P_i formulas used above:")
    for i in range(b + 1):
        print(f"P_{i} = {build_p(i, b)}")


if __name__ == "__main__":
    main()
