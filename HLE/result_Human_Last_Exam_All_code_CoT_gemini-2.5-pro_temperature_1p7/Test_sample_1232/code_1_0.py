import numpy as np

def check_log_supermodularity(q, beta):
    """
    Checks the log-supermodularity condition for the Potts model pair interaction.
    This condition is necessary for the FKG inequality to hold on a K_2 graph.
    The condition is: W(i,j) * W(i+1,j+1) >= W(i,j+1) * W(i+1,j) for all i, j.
    """
    if beta == 0:
        w = 1.0
    else:
        w = np.exp(2 * beta)

    # We only need to find one violation
    for i in range(1, q):       # i from 1 to q-1
        for j in range(1, q):   # j from 1 to q-1
            # W_ij = w if i==j else 1
            w_i_j = w if i == j else 1
            w_i1_j1 = w if (i + 1) == (j + 1) else 1
            w_i_j1 = w if i == (j + 1) else 1
            w_i1_j = w if (i + 1) == j else 1

            lhs = w_i_j * w_i1_j1
            rhs = w_i_j1 * w_i1_j
            
            if lhs < rhs:
                print(f"Violation found for q={q}, beta={beta:.2f}:")
                print(f"Condition: W({i},{j}) * W({i+1},{j+1}) >= W({i},{j+1}) * W({i+1},{j})")
                print(f"Values: {w_i_j} * {w_i1_j1} >= {w_i_j1} * {w_i1_j}")
                print(f"Equation: {lhs} >= {rhs} -> This is FALSE.")
                return False
                
    return True

# Main logic
# We are looking for the largest d. Let's test d=1.
# A graph with max_degree=1 is K_2. Let's test the condition for K_2.
# The condition must hold for ALL q >= 2.

# Test q=2
q2_holds = check_log_supermodularity(q=2, beta=1.0)
print(f"For q=2, the condition holds: {q2_holds}")
print("-" * 20)

# Test q=3
q3_holds = check_log_supermodularity(q=3, beta=1.0)
print(f"For q=3, the condition holds: {q3_holds}")
print("-" * 20)

print("The analysis shows:")
print("For q=2, the log-supermodularity holds, so the FKG inequality holds for any graph, and thus any max degree d.")
print("For q=3, the log-supermodularity fails for a single edge (a graph with max degree d=1).")
print("The question asks for the largest d where the property holds for ALL q >= 2.")
print("Since the property fails for d=1 with q=3, d cannot be 1 or higher.")
print("The only remaining possibility is d=0.")
print("A connected graph with max degree d=0 is a single vertex, for which the property holds trivially.")
print("The largest d is therefore 0.")
print("\nFinal calculation demonstrating the failure from the text above:")
# With q=3, choose i=1, c=2, b=3, d=2.
a, c = 1, 2 # spins at vertex u
b, d = 3, 2 # spins at vertex v
lhs = (max(a, c) == max(b, d)) + (min(a, c) == min(b, d))
rhs = (a == b) + (c == d)
print(f"I(max({a},{c}) = max({b},{d})) + I(min({a},{c}) = min({b},{d})) >= I({a}={b}) + I({c}={d})")
print(f"I({max(a,c)} = {max(b,d)}) + I({min(a,c)} = {min(b,d)}) >= {int(a==b)} + {int(c==d)}")
print(f"{int(max(a, c) == max(b, d))} + {int(min(a, c) == min(b, d))} >= {rhs}")
print(f"{lhs} >= {rhs}")
print(f"This is {'True' if lhs >= rhs else 'False'}.")
