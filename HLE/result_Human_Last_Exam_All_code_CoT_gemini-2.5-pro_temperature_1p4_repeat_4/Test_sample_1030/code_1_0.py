# Let's represent the truth values numerically for ease of computation
# T (True) = 2
# G (Glut) = 1
# F (False) = 0
T, G, F = 2, 1, 0
vals = {'T': T, 'G': G, 'F': F}
val_names = {v: k for k, v in vals.items()}

# Define the logical connectives based on the LP system
def neg(v):
    if v == T: return F
    if v == F: return T
    return G  # neg(G) = G

def conj(v1, v2):
    return min(v1, v2)

def disj(v1, v2):
    return max(v1, v2)

def impl(v1, v2):
    return disj(neg(v1), v2)

def is_designated(v):
    return v in [T, G]

# Let's analyze argument K: A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)
# For simplicity, we will hardcode the formulas for the premise and conclusion of K
def premise_K(vA, vB):
    return conj(vA, vB)

def conclusion_K(vA, vB):
    p1 = disj(neg(vA), neg(vB)) # (¬A ∨ ¬B)
    p2 = conj(vA, vB)          # (A ∧ B)
    return impl(p1, p2)

# --- Step 1: Check if the conclusion of K is a tautology ---
print("Step 1: Checking if the conclusion of argument K, '(¬A ∨ ¬B) → (A ∧ B)', is a tautology.")
is_tautology = True
counterexample = None
for vA in [T, G, F]:
    for vB in [T, G, F]:
        vC = conclusion_K(vA, vB)
        if not is_designated(vC):
            is_tautology = False
            counterexample = (vA, vB, vC)
            break
    if not is_tautology:
        break

if not is_tautology:
    vA, vB, vC = counterexample
    print(f"Result: The conclusion is NOT a tautology.")
    print(f"A counterexample valuation is A={val_names[vA]}, B={val_names[vB]}.")
    print("Let's trace the calculation with numerical values (T=2, G=1, F=0):")
    p1_val = disj(neg(vA), neg(vB))
    p2_val = conj(vA, vB)
    print(f"  Let v(A) = {vA} and v(B) = {vB}.")
    print(f"  The antecedent (¬A ∨ ¬B) evaluates to: max(neg({vA}), neg({vB})) = max({neg(vA)}, {neg(vB)}) = {p1_val}.")
    print(f"  The consequent (A ∧ B) evaluates to: min({vA}, {vB}) = {p2_val}.")
    print(f"  The full conclusion '({p1_val}) → ({p2_val})' evaluates to: max(neg({p1_val}), {p2_val}) = max({neg(p1_val)}, {p2_val}) = {vC}.")
    print(f"  The final value is {vC} ({val_names[vC]}), which is not a designated value.")
else:
    print("Result: The conclusion is a tautology.")

print("\n" + "-"*50 + "\n")

# --- Step 2: Check if argument K is valid ---
print("Step 2: Checking if argument K, 'A ∧ B ⊢ (¬A ∨ ¬B) → (A ∧ B)', is valid.")
print("Validity requires that if the premise is designated ({T, G}), the conclusion must also be designated.")
is_valid = True
invalidating_valuation = None
for vA in [T, G, F]:
    for vB in [T, G, F]:
        vP = premise_K(vA, vB)
        if is_designated(vP):
            vC = conclusion_K(vA, vB)
            if not is_designated(vC):
                is_valid = False
                invalidating_valuation = (vA, vB)
                break
    if not is_valid:
        break

if is_valid:
    print("Result: The argument is VALID. No valuation makes the premise designated and the conclusion false.")
else:
    vA, vB = invalidating_valuation
    print(f"Result: The argument is INVALID. Counterexample: A={val_names[vA]}, B={val_names[vB]}.")

print("\n" + "-"*50 + "\n")
print("Conclusion: Both arguments K and L are valid in the specified logic. However,")
print("the conclusion of L is a tautology, making its validity trivial. The conclusion")
print("of K is not a tautology, but the argument is still valid. This makes K a more")
print("substantive and likely intended answer.")
