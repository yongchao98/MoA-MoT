import collections

# --- Define the problem instance ---
# A multiset of natural numbers W
W = collections.Counter([1, 2, 3])
# The number of partitions
m = 2
# The target sum for each partition
b = 3

# --- Sanity check ---
if sum(k * v for k, v in W.items()) != m * b:
    print("Error: The sum of elements in W does not equal m * b.")
    exit()

# --- Define the base formulas (as strings) ---
# Using '1' for the constant 1 and 'bot' for the constant ⊥
# Using '*' for the connective ⊗ (tensor)
# Using '->' for the connective ⊸ (lollipop/linear implication)
A = "(1 -> bot)"
B = f"({A} -> {A})"

# --- Define the function f(w) ---
def generate_f(w_val):
    """Generates the formula string for f(w).
    f(w) = A -> (A -> ... (A -> B)...)  (w times)
    """
    res = B
    for _ in range(w_val):
        res = f"({A} -> {res})"
    return res

# --- Define the formula C(W, m, b) ---
def generate_C(W_set, m_val, b_val):
    """Generates the formula string for C(W, m, b)."""
    n = sum(W_set.values()) # Total number of items in the multiset W

    # A_b = A * A * ... * A (b times)
    A_b = " * ".join([A] * b_val)
    if b_val > 1:
        A_b = f"({A_b})"

    # A_b_m = A_b * A_b * ... * A_b (m times)
    A_b_m = " * ".join([A_b] * m_val)
    if m_val > 1:
        A_b_m = f"({A_b_m})"

    # B_n = B * B * ... * B (n times)
    B_n = " * ".join([B] * n)
    if n > 1:
        B_n = f"({B_n})"
    
    # consumer = B_n -> bot
    consumer = f"({B_n} -> bot)"

    # C = (A_b_m * consumer) -> bot
    C = f"(({A_b_m} * {consumer}) -> bot)"
    return C

# --- Print the results ---
print("--------------------------------------------------")
print("An encoding of the Equipartitioning Problem in Linear Logic")
print("--------------------------------------------------")
print(f"Problem Instance: W = {list(W.elements())}, m = {m}, b = {b}\n")
print("Base Formulas:")
print(f"  A = {A}")
print(f"  B = {B}\n")

print("Formulas f(w) for each w in W:")
for w in sorted(W.elements()):
    # We print each f(w) even for repeated numbers
    print(f"  f({w}) = {generate_f(w)}")

print("\nFormula C(W, m, b):")
C_formula = generate_C(W, m, b)
print(f"  C = {C_formula}")
print("\n--------------------------------------------------")
print("The sequent {f(w) | w in W} |- C is derivable iff EP(W, m, b) is true.")
print("--------------------------------------------------")
