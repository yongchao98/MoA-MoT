import math

def get_hamiltonicity_threshold():
    """
    This function provides the d-threshold for Hamiltonicity based on the problem description.

    The problem considers a graph H_n with minimum degree d = n/2 - eta, where
    1/2 <= eta <= n/64. We need to find the threshold probability p such that
    H_n U G(n, p) is Hamiltonian asymptotically almost surely.

    The threshold is determined by the 'worst-case' choice of H_n. For the given
    range of d, H_n can be constructed as a graph with a sparse cut, such as
    two disjoint cliques (for eta >= 1) or a complete bipartite graph (for eta=1/2).
    To make such graphs Hamiltonian, the random edges from G(n,p) must bridge
    this sparse cut.

    The analysis of these worst-case scenarios shows that a necessary condition
    (e.g., ensuring minimum degree 1 in the connecting random bipartite graph,
    or connecting a partition in the bipartite-like graph) requires a probability
    p of the order Omega(log(n)/n).

    An upper bound is given by the fact that G(n, p) itself is Hamiltonian for
    p = O(log(n)/n).

    Combining the lower and upper bounds, the d-threshold is Theta(log(n)/n) for the
    entire specified range of eta.
    """

    # The threshold p(n) is asymptotically equivalent to log(n)/n.
    # The problem asks to output the final equation.
    # The final equation is symbolic.
    
    equation_lhs = "p(n)"
    asymptotic_symbol = "="
    asymptotic_class = "Theta"
    equation_rhs_numerator_log = "log"
    equation_rhs_numerator_arg = "n"
    equation_rhs_denominator = "n"
    
    print("The d-threshold is:", end=" ")
    # The instruction "output each number in the final equation" is interpreted
    # as printing the components of the symbolic formula, as there are no explicit numbers.
    print(f"{equation_lhs} {asymptotic_symbol} {asymptotic_class}(({equation_rhs_numerator_log}({equation_rhs_numerator_arg}))/{equation_rhs_denominator})")

get_hamiltonicity_threshold()