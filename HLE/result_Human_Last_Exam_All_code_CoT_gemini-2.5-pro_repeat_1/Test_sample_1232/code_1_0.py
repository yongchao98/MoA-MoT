import itertools

def main():
    """
    Solves the problem by programmatically checking the conditions for the FKG inequality
    for the Potts model.
    """

    def I(condition):
        """Indicator function."""
        return 1 if condition else 0

    def check_local_fkg_for_potts(q):
        """
        Checks the local supermodularity condition for the q-state Potts model energy function.
        A violation for a given q means the FKG inequality fails for any graph with at least
        one edge, which has a maximum degree of at least 1.
        """
        print(f"--- Checking for q = {q} states ---")
        states = range(1, q + 1)
        # Iterate over all spin assignments (a,b,c,d) for (xi(u), xi(v), eta(u), eta(v)) on an edge (u,v).
        for a, b, c, d in itertools.product(states, repeat=4):
            E_xi_plus_eta = I(a == b) + I(c == d)
            E_join_plus_meet = I(max(a, c) == max(b, d)) + I(min(a, c) == min(b, d))

            if E_join_plus_meet < E_xi_plus_eta:
                print(f"Violation found for q = {q}!")
                print(f"Let spins on an edge be xi=(xi(u),xi(v)) and eta=(eta(u),eta(v)).")
                print(f"Counterexample: xi=({a}, {b}), eta=({c}, {d})")
                print(f"The supermodularity inequality is: E(xi v eta) + E(xi ^ eta) >= E(xi) + E(eta)")
                print(f"In this case: I(max({a},{c})==max({b},{d})) + I(min({a},{c})==min({b},{d})) >= I({a}=={b}) + I({c}=={d})")
                print(f"Plugging in values: {I(max(a, c) == max(b, d))} + {I(min(a, c) == min(b, d))} >= {I(a == b)} + {I(c == d)}")
                print(f"Resulting in: {E_join_plus_meet} >= {E_xi_plus_eta}, which is FALSE.")
                return False

        print(f"No violation found for q = {q}. The local condition holds.")
        return True

    print("Step 1: Analyze the condition for the positive correlations property (FKG inequality).")
    print("The property holds if the energy function is supermodular. For the Potts model, this")
    print("boils down to a local condition on each edge of the graph.\n")

    print("Step 2: Test this local condition for different numbers of states q.")
    # Test for q=2
    q2_holds = check_local_fkg_for_potts(2)
    print("\n")

    # Test for q=3
    q3_holds = check_local_fkg_for_potts(3)
    print("\n")

    print("Step 3: Draw conclusions from the test results.")
    print(f"For q=2 (Ising model), the condition holds? {q2_holds}. This is a known result. The property holds for any graph.")
    print(f"For q=3, the condition holds? {q3_holds}. A counterexample was found.")
    print("The question requires the property to hold for ANY q >= 2.")
    print("Our counterexample for q=3 occurred on a single edge. A graph consisting of just")
    print("one edge (like K_2) is connected and has a maximum degree d=1.")
    print("This means for d=1, we have found a graph (K_2) and a q (q=3) for which the property fails.")
    print("Therefore, the statement in the question is FALSE for d=1.")
    print("If the statement is false for d=1, it must also be false for any d > 1, because")
    print("any class of graphs with max_degree <= d (for d>1) contains the graph K_2.\n")

    print("Step 4: Analyze the only remaining possibility, d=0.")
    print("A finite, connected graph G with maximum degree d=0 must consist of a single vertex and no edges.")
    print("For such a graph, the sum over edges in the Potts model exponent is empty (zero).")
    print("The Gibbs measure is pi(xi) = 1/Z * exp(0), which is the uniform distribution on {1, ..., q}.")
    print("The FKG inequality for a uniform distribution on an ordered set is a famous result known as")
    print("Chebyshev's Sum Inequality, which states that E[fg] >= E[f]E[g] for any two increasing functions f and g.")
    print("Let's demonstrate with an example for q=3:")
    f = [1, 2, 5] # An increasing function on {1,2,3}
    g = [3, 6, 7] # Another increasing function on {1,2,3}
    n = 3
    E_f = sum(f) / n
    E_g = sum(g) / n
    E_fg = sum(f[i]*g[i] for i in range(n)) / n
    print(f"  Let f = {f}, g = {g}. E[f] = {E_f:.2f}, E[g] = {E_g:.2f}, E[fg] = {E_fg:.2f}")
    print(f"  The inequality E[fg] >= E[f]*E[g] becomes {E_fg:.2f} >= {E_f*E_g:.2f}, which is {E_fg >= E_f*E_g}.")
    print("So, the property holds for d=0.\n")

    print("Step 5: Final Conclusion.")
    print("The property holds for d=0.")
    print("The property fails for d=1 and any d > 1.")
    final_answer = 0
    print(f"The largest integer d for which the statement is true is {final_answer}.")


if __name__ == "__main__":
    main()