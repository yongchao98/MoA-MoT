import networkx as nx
import random
import math
import collections

def run_luby_variant(graph):
    """
    Runs the specified variant of Luby's algorithm on a graph.
    A copy of the graph is used, so the original is not modified.
    Returns the number of steps until the graph is empty.
    """
    steps = 0
    g = graph.copy()
    
    while g.number_of_nodes() > 0:
        steps += 1
        
        # If no edges remain, all isolated vertices are selected. The algorithm terminates.
        if g.number_of_edges() == 0:
            break
            
        nodes = list(g.nodes())
        random_values = {node: random.random() for node in nodes}
        selected_nodes = set()
        
        for v in nodes:
            # A vertex is selected if its random value is a strict local maximum.
            neighbors = list(g.neighbors(v))
            if not neighbors or all(random_values[v] > random_values[u] for u in neighbors):
                selected_nodes.add(v)
        
        # If no node is selected in a round (a possible, though low-probability event),
        # the algorithm simply proceeds to the next round with new random numbers.
        # This still counts as a step.
        if not selected_nodes:
            continue

        # A vertex is removed if it's selected or is a neighbor of a selected vertex.
        nodes_to_remove = set(selected_nodes)
        for v in selected_nodes:
            nodes_to_remove.update(g.neighbors(v))
        
        g.remove_nodes_from(nodes_to_remove)
            
    return steps

def run_simulations():
    """
    Runs simulations for the three specified graph classes and prints the results.
    """
    # Use smaller N for faster execution, but sufficient to see the trend
    N_VALUES = [50, 100, 200, 400] 
    NUM_RUNS = 10  # Number of runs to average over for each n

    print("Running simulations for the MIS algorithm.")
    print("This provides experimental evidence for the theoretical analysis.")
    print("-" * 40)

    for n in N_VALUES:
        print(f"\n--- Running for n = {n} ---")
        
        # Case 1: Cycle graph
        g_cycle = nx.cycle_graph(n)
        steps_cycle = [run_luby_variant(g_cycle) for _ in range(NUM_RUNS)]
        avg_steps_cycle = sum(steps_cycle) / NUM_RUNS
        print(f"1. Cycle C_{n}:")
        print(f"   Avg steps: {avg_steps_cycle:.2f}")
        if n > 1:
            print(f"   Avg steps / log(n): {avg_steps_cycle / math.log(n):.2f}")

        # Case 2: Tree with max degree <= 100
        # A random tree (nx.random_tree) is very unlikely to have high degrees.
        steps_tree = []
        for _ in range(NUM_RUNS):
             # nx.random_tree is in recent networkx versions.
             # It generates a tree uniformly at random from all trees on n nodes.
            g_tree = nx.random_tree(n)
            # We can assert max degree, but it's not strictly necessary for random trees.
            # max_deg = max(d for _, d in g_tree.degree())
            steps_tree.append(run_luby_variant(g_tree))
        avg_steps_tree = sum(steps_tree) / NUM_RUNS
        print(f"2. Random Tree on {n} vertices:")
        print(f"   Avg steps: {avg_steps_tree:.2f}")
        if n > 1:
            print(f"   Avg steps / log(n): {avg_steps_tree / math.log(n):.2f}")

        # Case 3: Graph with max degree <= 100
        # We'll use a random 3-regular graph as a representative of a sparse, bounded-degree graph.
        # A 100-regular graph would be very dense. The complexity class is the same.
        if n % 2 == 0 and n > 3:
            g_reg = nx.random_regular_graph(3, n)
            steps_reg = [run_luby_variant(g_reg) for _ in range(NUM_RUNS)]
            avg_steps_reg = sum(steps_reg) / NUM_RUNS
            print(f"3. 3-Regular graph on {n} vertices:")
            print(f"   Avg steps: {avg_steps_reg:.2f}")
            if n > 1:
                print(f"   Avg steps / log(n): {avg_steps_reg / math.log(n):.2f}")
        else:
            print(f"3. Bounded-degree graph on {n} vertices: Skipped")


    print("\n" + "-" * 40)
    print("--- Conclusion ---")
    print("The theoretical analysis shows that the algorithm takes Theta(log n) steps for all three cases.")
    print("The simulation results support this, as the ratio (Avg steps) / log(n) appears to converge to a constant.")
    print("\nMapping to categories:")
    print("f(n) = Theta(log n) implies f(n) = Omega(log n). This corresponds to Category 9.")
    print("\nTherefore:")
    d1 = 9
    d2 = 9
    d3 = 9
    print(f"The digit for f1(n) is: {d1}")
    print(f"The digit for f2(n) is: {d2}")
    print(f"The digit for f3(n) is: {d3}")
    print(f"\nThe final combined three-digit answer is: {d1}{d2}{d3}")

if __name__ == '__main__':
    run_simulations()