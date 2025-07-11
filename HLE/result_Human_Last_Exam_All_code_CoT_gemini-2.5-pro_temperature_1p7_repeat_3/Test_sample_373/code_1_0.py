import sys

def solve_barbell_gossip_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a
    barbell graph with 10 nodes during randomized uniform gossiping.
    """
    # Total number of nodes in the barbell graph
    total_nodes = 10

    # The barbell graph consists of two K_5 cliques.
    # So the size of each clique is 5.
    clique_size = 5

    # The nodes connected by the bottleneck edge are u and v.
    # The degree of such a node is the number of other nodes in its clique (clique_size - 1)
    # plus the connection to the other bottleneck node (+1).
    degree_bottleneck_node = (clique_size - 1) + 1

    # Probability of picking a specific node (e.g., u or v)
    # is 1 / total_nodes.
    prob_pick_node = 1 / total_nodes

    # Given a bottleneck node is picked, the probability of it choosing its
    # partner across the bridge is 1 / its degree.
    prob_pick_neighbor = 1 / degree_bottleneck_node

    # The total probability is the sum of two scenarios:
    # 1. Node u is picked and chooses v.
    # 2. Node v is picked and chooses u.
    # P(total) = P(pick u) * P(u chooses v) + P(pick v) * P(v chooses u)
    total_probability = (prob_pick_node * prob_pick_neighbor) + (prob_pick_node * prob_pick_neighbor)

    # Print the equation with all the numbers, as requested.
    # Redirecting output to stdout for the final answer
    original_stdout = sys.stdout
    sys.stdout = sys.__stdout__
    
    print(f"The probability of sampling the bottleneck edge is the sum of two scenarios:")
    print(f"1. Node 'u' is picked and it chooses node 'v': (1/{total_nodes}) * (1/{degree_bottleneck_node})")
    print(f"2. Node 'v' is picked and it chooses node 'u': (1/{total_nodes}) * (1/{degree_bottleneck_node})")
    print(f"\nTotal Probability = (1/{total_nodes} * 1/{degree_bottleneck_node}) + (1/{total_nodes} * 1/{degree_bottleneck_node})")
    print(f"Total Probability = {prob_pick_node * prob_pick_neighbor} + {prob_pick_node * prob_pick_neighbor} = {total_probability}")

    sys.stdout = original_stdout # Restore stdout if necessary, though it might not be in this simple script.

solve_barbell_gossip_probability()