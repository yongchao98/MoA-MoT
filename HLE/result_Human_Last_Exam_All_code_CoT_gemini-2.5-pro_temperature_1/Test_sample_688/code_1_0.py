import networkx as nx
from math import factorial
from itertools import combinations
from fractions import Fraction

def calculate_cn(n):
    """
    Calculates the prefactor c_n for a given n.

    The formula is:
    c_n = - (n-1)/n! * SUM_{E | K_n-E is biconnected} (-1)^|E|
    """
    if n < 2:
        return 0

    # The sum part of the formula
    the_sum = 0

    # Handle n=2 as a special case, as networkx.is_biconnected is for n>=3.
    # For n=2, K_2 is the only biconnected graph. It is K_2-E where E is the empty set.
    # The sum is (-1)^0 = 1.
    if n == 2:
        the_sum = 1
    else:
        kn = nx.complete_graph(n)
        edges = list(kn.edges())
        num_edges_kn = len(edges)
        
        # Iterate over all possible subsets of edges to remove
        for i in range(num_edges_kn + 1):
            for E in combinations(edges, i):
                # Create the subgraph H = K_n - E
                H = kn.copy()
                H.remove_edges_from(E)
                
                # Check if the resulting graph is biconnected.
                # nx.is_biconnected checks for connectivity internally.
                if nx.is_biconnected(H):
                    the_sum += (-1)**len(E)
    
    # Calculate the prefactor
    prefactor_num = -(n - 1)
    prefactor_den = factorial(n)
    
    # The result is a fraction
    result = Fraction(prefactor_num * the_sum, prefactor_den)
    return result

def main():
    """
    Calculates and prints c_n for n from 2 to 6.
    """
    print("Calculating the prefactor c_n for n = 2 to 6:")
    max_n = 6
    for n in range(2, max_n + 1):
        cn = calculate_cn(n)
        print(f"For n={n}, the prefactor is c_{n} = {cn.numerator}/{cn.denominator}")

if __name__ == "__main__":
    main()
