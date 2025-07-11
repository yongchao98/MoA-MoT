import collections

class Node:
    """A simple class to represent a node in the parse tree."""
    def __init__(self, value, layer=1, parent=None):
        self.value = value
        self.parent = parent
        self.children = []
        self.layer = layer

    def add_child(self, value):
        """Adds a child node and returns it."""
        child_node = Node(value, layer=self.layer + 1, parent=self)
        self.children.append(child_node)
        return child_node

    def __repr__(self):
        return f"Node(value='{self.value}', layer={self.layer})"

def build_parse_tree():
    """
    Builds the parse tree for the expression: y + (x + 4) * 5
    """
    # Layer 1
    root = Node("<expression>", layer=1)
    
    # Layer 2: <expression> -> <expression> + <term>
    e1 = root.add_child("<expression>")
    root.add_child("+")
    t1 = root.add_child("<term>")

    # Layer 3: Derivations for 'y' and '(x + 4) * 5'
    t2 = e1.add_child("<term>")                 # for y
    t3 = t1.add_child("<term>")                 # for (x + 4)
    t1.add_child("*")
    f1 = t1.add_child("<factor>")               # for 5

    # Layer 4
    f2 = t2.add_child("<factor>")               # for y
    f3 = t3.add_child("<factor>")               # for (x + 4)
    num1 = f1.add_child("number")               # for 5
    
    # Layer 5
    n1 = f2.add_child("name")                   # for y
    f3.add_child("(")
    e2 = f3.add_child("<expression>")           # for x + 4
    f3.add_child(")")
    num1.add_child("5")

    # Layer 6
    n1.add_child("y")
    e3 = e2.add_child("<expression>")           # for x
    e2.add_child("+")
    t4 = e2.add_child("<term>")                 # for 4

    # Layer 7
    t5 = e3.add_child("<term>")                 # for x
    f4 = t4.add_child("<factor>")               # for 4

    # Layer 8
    f5 = t5.add_child("<factor>")               # for x
    num2 = f4.add_child("number")               # for 4

    # Layer 9
    n2 = f5.add_child("name")                   # for x
    num2.add_child("4")

    # Layer 10
    n2.add_child("x")

    return root

def get_all_nodes(root):
    """Performs a breadth-first traversal to get all nodes."""
    q = collections.deque([root])
    nodes = []
    while q:
        node = q.popleft()
        nodes.append(node)
        for child in node.children:
            q.append(child)
    return nodes

def get_layers(all_nodes):
    """Groups all nodes by their layer number."""
    layers = collections.defaultdict(list)
    for node in all_nodes:
        layers[node.layer].append(node)
    return dict(sorted(layers.items()))

# --- Check Functions for Statements A-H ---

def check_A(all_nodes):
    # A. There is at least one <expression> which has a parent that is also an <expression> node.
    for node in all_nodes:
        if node.value == "<expression>" and node.parent and node.parent.value == "<expression>":
            return True
    return False

def check_B(layers):
    # B. The deepest number node is in the second to last layer of the tree.
    max_layer = max(layers.keys())
    deepest_number_layer = 0
    for node in get_all_nodes(list(layers[1])[0]):
        if node.value == "number":
            deepest_number_layer = max(deepest_number_layer, node.layer)
    return deepest_number_layer == (max_layer - 1)

def check_C(all_nodes):
    # C. There is a name node that appears in a layer which is between ... two layers ... that each contain a number node.
    number_layers = {node.layer for node in all_nodes if node.value == "number"}
    name_layers = {node.layer for node in all_nodes if node.value == "name"}
    for nl in name_layers:
        has_smaller = any(num_l < nl for num_l in number_layers)
        has_larger = any(num_l > nl for num_l in number_layers)
        if has_smaller and has_larger:
            return True
    return False

def check_D(layers):
    # D. The deepest layer contains a name with a <factor> as a parent.
    # Note: A 'name' node is a non-terminal. The deepest layer contains terminals.
    deepest_layer_num = max(layers.keys())
    for node in layers[deepest_layer_num]:
        # The premise is "The deepest layer contains a name [node]".
        if node.value == "name":
            # If a 'name' node is found, check its parent.
            if node.parent and node.parent.value == "<factor>":
                return True
    return False # No 'name' node found in the deepest layer.

def check_E(layers):
    # E. There is a layer that only has <factor> nodes, one operator, and one <term> node.
    operators = {'+', '-', '*', '/'}
    for layer_num, nodes_in_layer in layers.items():
        values = [n.value for n in nodes_in_layer]
        # Count node types in the current layer
        factor_count = values.count("<factor>")
        term_count = values.count("<term>")
        op_count = sum(1 for v in values if v in operators)
        
        # Check if the layer contains ONLY these types of nodes
        if (factor_count + term_count + op_count) != len(values):
            continue
        
        if factor_count >= 1 and term_count == 1 and op_count == 1:
            return True
    return False

def check_F(layers):
    # F. The node in the deepest layer has a parent which in turn has a <term> as a parent.
    deepest_layer_num = max(layers.keys())
    for node in layers[deepest_layer_num]:
        # Check for parent and grandparent
        if node.parent and node.parent.parent:
            if node.parent.parent.value == "<term>":
                return True
    return False

def check_G(layers):
    # G. There are 9 layers in the parse tree...
    return len(layers) == 9

def check_H(layers):
    # H. The fourth layer contains only <factor> nodes except for a single number node.
    if 4 not in layers:
        return False
    layer4_values = [n.value for n in layers[4]]
    number_count = layer4_values.count("number")
    factor_count = layer4_values.count("<factor>")
    return number_count == 1 and (number_count + factor_count) == len(layer4_values)

def main():
    # Build and analyze the tree
    tree_root = build_parse_tree()
    all_nodes = get_all_nodes(tree_root)
    layers = get_layers(all_nodes)
    
    # Evaluate all statements
    results = {
        'A': check_A(all_nodes),
        'B': check_B(layers),
        'C': check_C(all_nodes),
        'D': check_D(layers),
        'E': check_E(layers),
        'F': check_F(layers),
        'G': check_G(layers),
        'H': check_H(layers),
    }

    print("Analysis of each statement:")
    false_statements = []
    for statement, is_true in results.items():
        print(f"Statement {statement} is {is_true}.")
        if not is_true:
            false_statements.append(statement)

    print("\nConclusion:")
    if len(false_statements) > 1:
        print(f"Statements {', '.join(false_statements)} are false.")
        print("Since more than one statement is false, the correct option is I.")
        final_answer = 'I'
    elif len(false_statements) == 1:
        print(f"The only false statement is {false_statements[0]}.")
        final_answer = false_statements[0]
    else:
        print("All statements are true.")
        final_answer = 'None' # Should not happen based on the question

    print(f"\n<<<{final_answer}>>>")

if __name__ == "__main__":
    main()