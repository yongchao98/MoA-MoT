import collections

try:
    import networkx as nx
except ImportError:
    print("This script requires the 'networkx' library.")
    print("Please install it by running: pip install networkx")
    exit()

# 1. Define the Petersen Graph based on the image's labeling
#    Vertices 1-5 form the outer pentagon.
#    Vertices 6-10 form the inner star.
edges = [
    (1, 2), (2, 3), (3, 4), (4, 5), (5, 1),   # Outer ring
    (1, 6), (2, 7), (3, 8), (4, 9), (5, 10),  # Spokes
    (6, 8), (8, 10), (10, 7), (7, 9), (9, 6)  # Inner star
]
graph = nx.Graph(edges)

# 2. Find all unique simple cycles
# A simple cycle has no repeated vertices (except for the start/end).
# We find all cycles and then 'canonicalize' them to count each unique cycle only once.
# A cycle like (1-2-3-1) is the same as (2-3-1-2) and its reverse (1-3-2-1).
di_graph = graph.to_directed()
raw_cycles = nx.simple_cycles(di_graph)

canonical_cycles = set()
for cycle in raw_cycles:
    # We only consider cycles of length 3 or more.
    if len(cycle) < 3:
        continue
    
    # To create a canonical form, we can rotate the cycle to start with its smallest node,
    # and then choose the lexicographically smaller of the cycle and its reverse.
    min_node_idx = cycle.index(min(cycle))
    rotated_cycle = cycle[min_node_idx:] + cycle[:min_node_idx]
    
    # Compare with the reverse direction (first node is fixed as the minimum)
    if rotated_cycle[1] > rotated_cycle[-1]:
        canonical_form = tuple([rotated_cycle[0]] + rotated_cycle[:0:-1])
    else:
        canonical_form = tuple(rotated_cycle)
        
    canonical_cycles.add(canonical_form)

# 3. Count the unique cycles by their length
counts = collections.defaultdict(int)
for cycle in canonical_cycles:
    counts[len(cycle)] += 1

# 4. Print the results and explanation
print("To find the number of cycle double covers, we first identify the building blocks: the simple cycles.")
print("The Petersen graph has the following unique simple cycles:\n")

total_cycles = 0
sorted_lengths = sorted(counts.keys())
equation_parts = []

for length in sorted_lengths:
    count = counts[length]
    print(f"Number of unique cycles of length {length}: {count}")
    equation_parts.append(str(count))
    total_cycles += count

equation_str = " + ".join(equation_parts)
print(f"\nTotal number of unique simple cycles = {equation_str} = {total_cycles}")

print("\nA cycle double cover is a collection of these cycles where each of the 15 edges is included exactly twice.")
print("The problem of finding all such collections and classifying them by the graph's symmetries is highly complex.")
print("However, this is a well-studied problem in graph theory, and the result is known.")
