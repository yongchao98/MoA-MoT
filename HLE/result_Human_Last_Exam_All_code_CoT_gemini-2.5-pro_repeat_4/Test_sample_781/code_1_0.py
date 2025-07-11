import math

def solve():
    """
    Calculates the largest number n based on the Erdos-Gallai theorem for
    the maximum number of edges in a graph with a given matching number.
    """
    
    # v is the number of special points in the continuum X
    v = 5
    # k is the maximum size of a matching in the corresponding graph G.
    # The topological constraints imply that the matching number is at most 1.
    k = 1

    print(f"The problem reduces to finding the maximum number of edges in a graph with v={v} vertices and matching number at most k={k}.")
    print("We use the Erdos-Gallai Theorem for this calculation.")
    print(f"The formula is: max( C(2*k+1, 2), C(k, 2) + k*(v-k) )")
    print("-" * 20)

    # Calculate the first term of the max function
    # This corresponds to the number of edges in a complete graph K_{2k+1}
    arg1_val1 = 2 * k + 1
    term1 = math.comb(arg1_val1, 2)
    print(f"First term: C(2*{k}+1, 2) = C({arg1_val1}, 2) = {term1}")
    print("This corresponds to a K3 graph (a triangle), which has 3 edges and matching number 1.")

    # Calculate the second term of the max function
    # This corresponds to a star-like graph structure
    arg2_val1 = math.comb(k, 2)
    arg2_val2 = k * (v - k)
    term2 = arg2_val1 + arg2_val2
    print(f"Second term: C({k}, 2) + {k}*({v}-{k}) = {arg2_val1} + {arg2_val2} = {term2}")
    print("This corresponds to a star graph, where one vertex is connected to the other 4.")
    print("-" * 20)
    
    # The result is the maximum of the two terms
    result = max(term1, term2)
    print(f"The largest number n is the maximum of these two values.")
    print(f"n = max({term1}, {term2}) = {result}")

solve()