import itertools

def check_holley_for_potts(q):
    """
    Checks if the Holley criterion holds for the q-state Potts model
    with the standard linear order on the states {1, ..., q}.

    The criterion simplifies to a local condition for each edge:
    I{max(a,c)=max(b,d)} + I{min(a,c)=min(b,d)} >= I{a=b} + I{c=d}
    for all a,b,c,d in {1, ..., q}, where a,b are the spin values
    of two adjacent vertices for a configuration xi, and c,d are the
    spin values for a configuration eta.

    Args:
        q (int): The number of states in the Potts model.

    Returns:
        tuple: (bool, tuple or None)
               A tuple containing a boolean indicating if the criterion holds,
               and the first found counterexample if it does not.
    """
    states = range(1, q + 1)
    # itertools.product generates all q^4 combinations for (a, b, c, d)
    for a, b, c, d in itertools.product(states, repeat=4):
        # Left-hand side of the inequality
        lhs = (1 if max(a, c) == max(b, d) else 0) + \
              (1 if min(a, c) == min(b, d) else 0)

        # Right-hand side of the inequality
        rhs = (1 if a == b else 0) + \
              (1 if c == d else 0)

        if lhs < rhs:
            # Found a counterexample
            return False, (a, b, c, d)

    # No counterexample found
    return True, None

def main():
    """
    Main function to analyze the positive correlations property for the Potts model.
    """
    print("Analyzing the positive correlations property for the q-state Potts model.")
    print("The property holds if a local condition (from the Holley criterion) is met for all edges.")
    print("This condition depends on the number of states, q.")
    print("-" * 70)

    # Check for q=2 (Ising model equivalent)
    q2 = 2
    print(f"Checking for q = {q2}...")
    holds_q2, _ = check_holley_for_potts(q2)
    if holds_q2:
        print(f"Result: For q = {q2}, the condition holds for all spin combinations.")
        print("This implies the property holds for any graph G and any beta >= 0, when q = 2.")
    else:
        print(f"Result: For q = {q2}, a counterexample was found.")

    print("-" * 70)

    # Check for q=3
    q3 = 3
    print(f"Checking for q = {q3}...")
    holds_q3, counterexample_q3 = check_holley_for_potts(q3)
    if not holds_q3:
        print(f"Result: For q = {q3}, a counterexample was found.")
        a, b, c, d = counterexample_q3
        lhs = (1 if max(a, c) == max(b, d) else 0) + (1 if min(a, c) == min(b, d) else 0)
        rhs = (1 if a == b else 0) + (1 if c == d else 0)
        
        print(f"\nLet spins on an edge be (xi(x), xi(y)) = ({a}, {b}) and (eta(x), eta(y)) = ({c}, {d}).")
        print("The inequality to check is:")
        print(f"I{{max({a},{c})==max({b},{d})}} + I{{min({a},{c})==min({b},{d})}} >= I{{{a}=={b}}} + I{{{c}=={d}}}")
        print("Plugging in the numbers:")
        print(f"I{{{max(a,c)}=={max(b,d)}}} + I{{{min(a,c)}=={min(b,d)}}} >= {1 if a==b else 0} + {1 if c==d else 0}")
        print(f"{'1' if max(a,c)==max(b,d) else '0'} + {'1' if min(a,c)==min(b,d) else '0'} >= {rhs}")
        print(f"{lhs} >= {rhs}  <-- This is FALSE.")
        
        print("\nThis counterexample shows that for q >= 3, the property fails for any graph with at least one edge")
        print("if beta > 0, because the Holley criterion is violated.")

    print("-" * 70)
    print("The question asks for the largest maximum degree d such that the property holds for ANY q >= 2.")
    print("Since the property fails for q=3 on any graph with an edge (i.e., deg_max >= 1),")
    print("the statement cannot be true for d=1, 2, 3, etc.")
    print("We must therefore consider d=0. A connected graph with deg_max=0 is a single vertex.")
    print("For a single vertex graph, the property holds for all q and beta.")
    
    final_answer = 0
    print(f"\nConclusion: The largest integer d for which the statement is true is {final_answer}.")

if __name__ == "__main__":
    main()