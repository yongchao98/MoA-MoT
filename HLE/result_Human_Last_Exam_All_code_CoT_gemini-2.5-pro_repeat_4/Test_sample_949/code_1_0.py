import math
import sys

# Increase recursion limit for deep strings
sys.setrecursionlimit(2000)

class Node:
    """Represents a pair of parentheses as a node in a tree."""
    def __init__(self):
        self.children = []
        self.L = 0  # Length
        self.D = 0  # Depth

def parse_to_tree(s: str):
    """Parses a parenthesis string into a tree structure."""
    if not s:
        return None
    root = Node()
    stack = [root]
    # This simplified parser assumes a single root element wrapping the string
    # e.g., for "()()", it processes it as two children of an implicit root.
    # To handle multiple top-level pairs, we add a virtual root.
    s = f"({s})"
    
    node_map = {} # map node to its start index
    open_indices_stack = []

    # First pass to build tree structure
    root = Node()
    current_parent = root
    parent_stack = []

    i = 1 # Skip the virtual opening parenthesis
    while i < len(s) - 1:
        char = s[i]
        if char == '(':
            new_node = Node()
            current_parent.children.append(new_node)
            parent_stack.append(current_parent)
            current_parent = new_node
            i += 1
        elif char == ')':
            current_parent = parent_stack.pop()
            i += 1
    
    # Second pass (post-order traversal) to compute L and D
    all_pairs = []
    def compute_L_D(node: Node):
        if not node.children:
            node.D = 1
            node.L = 2
        else:
            max_child_D = 0
            sum_child_L = 0
            for child in node.children:
                compute_L_D(child)
                if child.D > max_child_D:
                    max_child_D = child.D
                sum_child_L += child.L
            node.D = 1 + max_child_D
            node.L = 2 + sum_child_L
        all_pairs.append({'L': node.L, 'D': node.D})

    # The virtual root's children are the actual top-level pairs
    for child_node in root.children:
        compute_L_D(child_node)

    return all_pairs

def flog(x):
    """Logarithm used in the problem, max(1, log(x))."""
    if x <= 1: return 1.0
    return max(1.0, math.log(x))

def evaluate_statements(pairs):
    """Calculates the sums for all six statements."""
    sums = {
        'L1': 0, 'D1': 0, 'L2': 0, 'D2': 0, 'L3': 0, 'D3': 0,
        'L4': 0, 'D4': 0, 'L5': 0, 'D5': 0, 'L6': 0, 'D6': 0,
    }
    for p in pairs:
        L, D = p['L'], p['D']
        logL, logD = flog(L), flog(D)
        
        sums['L1'] += logL
        sums['D1'] += logD
        
        sums['L2'] += flog(logL)
        sums['D2'] += flog(logD)

        sums['L3'] += logL**5
        sums['D3'] += logD**5

        if logL > 0:
            sums['L4'] += 2**math.sqrt(logL)
        if logD > 0:
            sums['D4'] += 2**math.sqrt(logD)
            
        sums['L5'] += L**0.1
        sums['D5'] += D**0.11
        
        sums['L6'] += L**0.25
        sums['D6'] += D**0.5
        
    ratios = {}
    for i in range(1, 7):
        # Avoid division by zero, though D is always >= 1
        ratios[i] = sums[f'L{i}'] / sums[f'D{i}'] if sums[f'D{i}'] != 0 else float('inf')
        
    return ratios

def gen_binary_tree_str(k):
    """Generates strings from the S_k = (S_{k-1}S_{k-1}) family."""
    if k == 0:
        return "()"
    s_prev = gen_binary_tree_str(k - 1)
    return f"({s_prev}{s_prev})"

def gen_comet_str(m, k):
    """Generates strings from the comet family."""
    head = "()" * k
    tail = f"({head})"
    for _ in range(m):
        tail = f"({tail})"
    return tail

def main():
    """
    Main function to run the analysis and print the final answer.
    The code demonstrates the diverging and converging nature of the ratios
    for different string families, providing computational evidence for the answers.
    """
    
    print("--- Analysis for TFFFTT ---")
    print("\nDemonstration for Statements 2, 3, 4 (False): Using 'Binary Tree' family S_k=(S_{k-1}S_{k-1})")
    print("Ratios should diverge.")
    print("k | Ratio 2 (loglog) | Ratio 3 (log^5) | Ratio 4 (2^sqrt(log))")
    print("-" * 65)
    for k in range(2, 8):
        s = gen_binary_tree_str(k)
        pairs = parse_to_tree(s)
        ratios = evaluate_statements(pairs)
        print(f"{k:1d} | {ratios[2]:16.2f} | {ratios[3]:15.2e} | {ratios[4]:22.2e}")

    print("\nDemonstration for Statement 5 (False): Using 'Comet' family C_{m,k}")
    print("Ratio should diverge. We set k=m^a where a is between 1.11 and 10/9.")
    print("m |     k     | Ratio 5 (L^0.1 / D^0.11)")
    print("-" * 45)
    a = 1.111 # This value is in [1.11, 10/9)
    for m in [5, 10, 15, 20]:
        k = int(m**a)
        s = gen_comet_str(m, k)
        pairs = parse_to_tree(s)
        ratios = evaluate_statements(pairs)
        print(f"{m:2d}| {k:9d} | {ratios[5]:25.4f}")


    print("\nDemonstration for Statements 1, 6 (True): Using 'Comet' family C_{m,k}")
    print("Ratios should converge/be bounded. Test with k=m^2.")
    print("m |      k     | Ratio 1 (log) | Ratio 6 (L^0.25 / D^0.5)")
    print("-" * 60)
    a = 2 
    for m in [5, 10, 15, 20]:
        k = int(m**a)
        s = gen_comet_str(m, k)
        pairs = parse_to_tree(s)
        ratios = evaluate_statements(pairs)
        print(f"{m:2d}| {k:10d} | {ratios[1]:13.4f} | {ratios[6]:25.4f}")

    # Based on theoretical analysis, the final answer is TFFFTT
    final_answer = "TFFFTT"
    print(f"\nFinal Answer String: {final_answer}")
    print(f"<<<{final_answer}>>>")

if __name__ == "__main__":
    main()
