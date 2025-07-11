import numpy as np

def heuristic(h1, h2):
    """
    Computes the heuristic [h1, h2, |h1-h2|, h1⊙h2] for numpy array inputs.
    """
    h1 = np.array(h1)
    h2 = np.array(h2)
    diff = np.abs(h1 - h2)
    hadamard = h1 * h2
    return np.concatenate([h1, h2, diff, hadamard])

# Define the mixed-dimension logical operators
def op_C_prime(p1, q2): return p1 and q2
def op_D_prime(p1, q2): return p1 or q2
def op_X_prime(p1, q2): return p1 != q2
def op_E_prime(p1, q2): return p1 == q2
def op_I_prime(p1, q2): return not p1 or q2 # p1 -> q2

operators_to_test = {
    "C'": op_C_prime,
    "D'": op_D_prime,
    "X'": op_X_prime,
    "E'": op_E_prime,
    "I'": op_I_prime
}

# We select four specific 2D input points (p1,q1) and (p2,q2)
# P1: (p1,q1,p2,q2) = (1,0,0,1)
# P2: (p1,q1,p2,q2) = (1,0,0,0)
# P3: (p1,q1,p2,q2) = (0,0,0,1)
# P4: (p1,q1,p2,q2) = (0,0,0,0)
h1_P1, h2_P1 = ([1, 0], [0, 1])
h1_P2, h2_P2 = ([1, 0], [0, 0])
h1_P3, h2_P3 = ([0, 0], [0, 1])
h1_P4, h2_P4 = ([0, 0], [0, 0])

# Compute the feature vectors for these points
feat_P1 = heuristic(h1_P1, h2_P1)
feat_P2 = heuristic(h1_P2, h2_P2)
feat_P3 = heuristic(h1_P3, h2_P3)
feat_P4 = heuristic(h1_P4, h2_P4)

print("--- Step 1: Find Linear Dependence in Feature Space ---")
# The four points are chosen such that their feature vectors are linearly dependent.
# Specifically, feat_P1 = feat_P2 + feat_P3 - feat_P4
dependency_check = np.array_equal(feat_P1, feat_P2 + feat_P3 - feat_P4)
print(f"Checking if feat(P1) = feat(P2) + feat(P3) - feat(P4): {dependency_check}")
print("This dependency implies that for any linear model score(h) = w•h + b,")
print("score(P1) = score(P2) + score(P3) - score(P4) + b.")
print("Since score(P4) = w•feat(P4) + b = w•0 + b = b, this simplifies to:")
print("score(P1) = score(P2) + score(P3)")
print("\n--- Step 2: Test Learnability of Each Operator ---")
print("An operator is NOT learnable if this equation leads to a contradiction.")
print("Contradiction: A positive number cannot be the sum of two negative numbers, and vice versa.")
print("-" * 50)

unlearnable_operators = []
for name, op_func in operators_to_test.items():
    print(f"Testing Operator: {name}")
    
    # Get the p1, q2 values for each point
    p1_P1, q2_P1 = h1_P1[0], h2_P1[1]
    p1_P2, q2_P2 = h1_P2[0], h2_P2[1]
    p1_P3, q2_P3 = h1_P3[0], h2_P3[1]
    
    # Calculate the target labels (1 for True, 0 for False)
    z1 = op_func(p1_P1, q2_P1)
    z2 = op_func(p1_P2, q2_P2)
    z3 = op_func(p1_P3, q2_P3)

    # Determine the required sign of the score for each point
    s1_sign = "positive (>0)" if z1 else "negative (<0)"
    s2_sign = "positive (>0)" if z2 else "negative (<0)"
    s3_sign = "positive (>0)" if z3 else "negative (<0)"
    
    # Check for contradiction
    contradiction = (z1 and not z2 and not z3) or (not z1 and z2 and z3)
    
    print(f"Labels: z(P1)={int(z1)}, z(P2)={int(z2)}, z(P3)={int(z3)}")
    print(f"Equation: score(P1) = score(P2) + score(P3)")
    print(f"Required signs: {s1_sign} = {s2_sign} + {s3_sign}")

    if contradiction:
        print("Result: CONTRADICTION. This operator is NOT learnable.")
        unlearnable_operators.append(name)
    else:
        print("Result: Consistent. This operator is potentially learnable.")
    print("-" * 50)

print("\n--- Conclusion ---")
print("All element-wise operators (X, C, D, E, I) are learnable.")
print("The mixed-dimension operators that are NOT learnable are:")
print(f"{{{', '.join(sorted(unlearnable_operators))}}}")
<<<K>>>