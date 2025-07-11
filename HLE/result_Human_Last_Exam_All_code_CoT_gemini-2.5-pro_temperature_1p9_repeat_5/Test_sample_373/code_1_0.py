import sys

# Suppress creation of __pycache__
sys.dont_write_bytecode = True

def solve():
    """
    Calculates the probability of sampling the bottleneck edge in a
    barbell graph with 10 nodes and uniform weights using randomized
    uniform gossiping.
    """

    # 1. Define the Graph
    # Total number of nodes in the barbell graph
    N = 10
    # A barbell graph consists of two equal-sized cliques connected by a bridge.
    # Size of each clique
    m = N / 2

    # Let the bottleneck edge connect node 'u' and node 'v'.

    # 2. Calculate Node Degrees
    # The degree of a node at the end of the bottleneck edge ('u' or 'v') is
    # the sum of its connections within its clique (m-1) and the bridge connection (1).
    degree_bottleneck_node = (m - 1) + 1

    # 3. & 4. Calculate the Probability
    # The probability of sampling an undirected edge {u, v} is the sum of probabilities
    # of sampling the directed edge (u, v) and the directed edge (v, u).
    # P({u,v}) = P(picking u) * P(picking v | u) + P(picking v) * P(picking u | v)
    # P({u,v}) = (1/N) * (1/deg(u)) + (1/N) * (1/deg(v))

    prob_u_then_v = (1/N) * (1/degree_bottleneck_node)
    prob_v_then_u = (1/N) * (1/degree_bottleneck_node)

    total_prob = prob_u_then_v + prob_v_then_u

    # 5. Print the results
    print("Let N be the total number of nodes and m be the size of each clique.")
    print(f"N = {N}, m = N/2 = {int(m)}\n")
    print("Let the bottleneck edge be {u, v}.")
    print("The degree of a node at the end of the bottleneck edge is (m-1) + 1.")
    print(f"deg(u) = deg(v) = ({int(m)}-1) + 1 = {int(degree_bottleneck_node)}\n")
    print("The probability of sampling the bottleneck edge {u,v} is:")
    print("P({u,v}) = (1/N) * (1/deg(u)) + (1/N) * (1/deg(v))")
    print(f"P({{u,v}}) = (1/{N}) * (1/{int(degree_bottleneck_node)}) + (1/{N}) * (1/{int(degree_bottleneck_node)})")
    print(f"P({{u,v}}) = {prob_u_then_v} + {prob_v_then_u}")
    print(f"P({{u,v}}) = {total_prob}")


solve()
<<<0.04>>>