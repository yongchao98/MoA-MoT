import numpy as np

def solve_d_zero_case():
    """
    This function analyzes the d=0 case for the problem.
    For d=0, the graph is a single vertex, and the Potts measure is uniform.
    The positive correlation property then reduces to Chebyshev's Sum Inequality.
    This code demonstrates that inequality for randomly generated increasing functions.
    """
    print("Analyzing the case where the maximum degree d = 0.")
    print("A connected graph with d=0 is a single vertex with no edges.")
    print("The Potts measure on such a graph is the uniform measure over the q states.")
    print("The positive correlation property must hold. We verify this using an example.")
    print("-" * 50)

    # Number of states in the Potts model
    q = 8

    # For the uniform measure on {1, ..., q}, an increasing function is
    # simply a non-decreasing sequence of q numbers.
    # Let's generate two random non-decreasing sequences to represent f and g.
    f_values = sorted(np.random.randint(1, 100, q))
    g_values = sorted(np.random.randint(1, 100, q))

    print(f"Let q = {q}")
    print(f"Let f be an increasing function with values f(1),...,f(q): {f_values}")
    print(f"Let g be another increasing function with values g(1),...,g(q): {g_values}")
    print("")

    # The positive correlation property is: E[fg] >= E[f]E[g]
    # For a uniform measure, E[h] = (1/q) * sum(h_i)
    E_f = np.mean(f_values)
    E_g = np.mean(g_values)
    
    fg_values = np.array(f_values) * np.array(g_values)
    E_fg = np.mean(fg_values)
    
    covariance = E_fg - (E_f * E_g)

    print("We check the inequality: E[fg] >= E[f] * E[g]")
    print(f"E[f] = sum(f) / q = {sum(f_values)} / {q} = {E_f}")
    print(f"E[g] = sum(g) / q = {sum(g_values)} / {q} = {E_g}")
    print(f"E[f] * E[g] = {E_f} * {E_g} = {E_f * E_g}")
    print("")
    print(f"E[fg] = sum(f*g) / q = {sum(fg_values)} / {q} = {E_fg}")
    print("")
    
    print(f"The inequality is {E_fg} >= {E_f * E_g}, which evaluates to {E_fg >= E_f * E_g}.")
    print(f"The covariance is E[fg] - E[f]E[g] = {covariance}, which is non-negative.")
    print("-" * 50)
    print("This property (Chebyshev's inequality) holds for any two non-decreasing sequences.")
    print("Thus, for d=0, positive correlation holds for any q.")
    print("\nFor any d >= 1, a graph with an edge exists, and a counterexample can be found for q >= 3.")
    print("Therefore, the largest d for which the statement is true is 0.")

solve_d_zero_case()