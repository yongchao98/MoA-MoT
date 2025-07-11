import itertools

def check_supermodularity_potts(q):
    """
    Checks if the supermodularity condition for the Potts Hamiltonian holds for a given q.
    The condition is: I(max(a,c)=max(b,d)) + I(min(a,c)=min(b,d)) >= I(a=b) + I(c=d)
    for all a,b,c,d in {1, ..., q}.
    """
    spins = range(1, q + 1)
    # itertools.product generates all possible combinations of a,b,c,d
    for a, b, c, d in itertools.product(spins, repeat=4):
        lhs = (max(a, c) == max(b, d)) + (min(a, c) == min(b, d))
        rhs = (a == b) + (c == d)
        if lhs < rhs:
            print(f"Supermodularity fails for q={q}.")
            print(f"Counterexample: a={a}, b={b}, c={c}, d={d}")
            print(f"LHS = I(max({a},{c})==max({b},{d})) + I(min({a},{c})==min({b},{d})) = {lhs}")
            print(f"RHS = I({a}=={b}) + I({c}=={d}) = {rhs}")
            print(f"Inequality {lhs} >= {rhs} is FALSE.")
            print("This failure occurs on a graph with a single edge, which has a maximum degree of 1.")
            return False
    print(f"Supermodularity holds for q={q}.")
    return True

print("Step 1: Checking the case q=2 (Ising Model)")
check_supermodularity_potts(2)
print("\n" + "="*50 + "\n")
print("Step 2: Checking the case q=3")
check_supermodularity_potts(3)
print("\nConclusion: The condition fails for q=3 on a graph with max_degree=1.")
print("Since the property must hold for ALL q>=2, any d>=1 is ruled out.")
print("For d=0, the graph is a single vertex, and the property holds (by Chebyshev's sum inequality).")
print("Therefore, the largest d is 0.")

# The final answer is an integer.
print("\nFinal Answer: 0")