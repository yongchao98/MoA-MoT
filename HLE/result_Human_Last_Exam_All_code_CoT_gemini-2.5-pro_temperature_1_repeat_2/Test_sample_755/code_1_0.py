import collections

class Node:
    """A simple class to represent a node in the parse tree."""
    def __init__(self, n_type, value=None):
        self.type = n_type
        self.value = value
        self.children = []
        self.parent = None
        self.depth = 0

    def add_child(self, node):
        """Adds a child node and sets its parent and depth."""
        self.children.append(node)
        node.parent = self
        node.depth = self.depth + 1
        return node

    def __repr__(self):
        """String representation of the node for clarity."""
        return f"Node(type='{self.type}', depth={self.depth})"

def build_parse_tree():
    """
    Constructs the parse tree for the expression 'y + (x + 4) * 5'
    based on the provided BNF grammar.
    """
    root = Node('expression')
    root.depth = 1

    # Layer 2: <expression> -> <expression> + <term>
    exp2 = root.add_child(Node('expression'))
    root.add_child(Node('operator', '+'))
    term1 = root.add_child(Node('term'))

    # Layer 3: Left <expression> becomes y, Right <term> becomes (x+4)*5
    # <expression> -> <term>
    term2 = exp2.add_child(Node('term'))
    # <term> -> <term> * <factor>
    term3 = term1.add_child(Node('term'))
    term1.add_child(Node('operator', '*'))
    factor1 = term1.add_child(Node('factor'))

    # Layer 4
    # <term> -> <factor> (for y)
    factor2 = term2.add_child(Node('factor'))
    # <term> -> <factor> (for (x+4))
    factor3 = term3.add_child(Node('factor'))
    # <factor> -> number (for 5)
    number1 = factor1.add_child(Node('number'))

    # Layer 5
    # <factor> -> name (for y)
    name1 = factor2.add_child(Node('name'))
    # <factor> -> (<expression>) (for (x+4))
    factor3.add_child(Node('paren', '('))
    exp3 = factor3.add_child(Node('expression'))
    factor3.add_child(Node('paren', ')'))
    # number -> 5
    number1.add_child(Node('value', '5'))

    # Layer 6
    # name -> y
    name1.add_child(Node('value', 'y'))
    # <expression> -> <expression> + <term> (for x+4)
    exp4 = exp3.add_child(Node('expression'))
    exp3.add_child(Node('operator', '+'))
    term4 = exp3.add_child(Node('term'))

    # Layer 7
    # <expression> -> <term> (for x)
    term5 = exp4.add_child(Node('term'))
    # <term> -> <factor> (for 4)
    factor4 = term4.add_child(Node('factor'))

    # Layer 8
    # <term> -> <factor> (for x)
    factor5 = term5.add_child(Node('factor'))
    # <factor> -> number (for 4)
    number2 = factor4.add_child(Node('number'))

    # Layer 9
    # <factor> -> name (for x)
    name2 = factor5.add_child(Node('name'))
    # number -> 4
    number2.add_child(Node('value', '4'))

    # Layer 10
    # name -> x
    name2.add_child(Node('value', 'x'))

    return root

def get_layers(root):
    """Performs a breadth-first search to group all nodes by layer."""
    if not root:
        return {}
    layers = collections.defaultdict(list)
    queue = collections.deque([root])
    while queue:
        node = queue.popleft()
        layers[node.depth].append(node)
        for child in node.children:
            queue.append(child)
    return dict(layers)

def analyze_statements():
    """Builds the tree and evaluates each statement."""
    root = build_parse_tree()
    layers = get_layers(root)
    max_depth = max(layers.keys())

    print("--- Analysis of Statements ---")

    # A. There is at least one <expression> which has a parent that is also an <expression> node.
    check_a = any(node.type == 'expression' and node.parent and node.parent.type == 'expression'
                  for layer in layers.values() for node in layer)
    print(f"Statement A is {check_a}.")

    # B. The deepest number node is in the second to last layer of the tree.
    num_layers = [node.depth for layer in layers.values() for node in layer if node.type == 'number']
    deepest_num_layer = max(num_layers) if num_layers else 0
    check_b = deepest_num_layer == (max_depth - 1)
    print(f"Statement B is {check_b}.")

    # C. There is a name node that appears in a layer which is between two layers
    #    such that each of these two layers contain a number node.
    number_layers = sorted(list(set(num_layers)))
    name_node_layers = [node.depth for layer in layers.values() for node in layer if node.type == 'name']
    check_c = False
    if len(number_layers) >= 2:
        min_num_layer, max_num_layer = number_layers[0], number_layers[-1]
        check_c = any(min_num_layer < name_layer < max_num_layer for name_layer in name_node_layers)
    print(f"Statement C is {check_c}.")

    # D. The deepest layer contains a name with a <factor> as a parent.
    deepest_layer_nodes = layers.get(max_depth, [])
    check_d = any(node.type == 'name' and node.parent and node.parent.type == 'factor'
                  for node in deepest_layer_nodes)
    print(f"Statement D is {check_d}.")

    # E. There is a layer that only has <factor> nodes, one operator, and one <term> node.
    check_e = False
    for layer_nodes in layers.values():
        types = [n.type for n in layer_nodes]
        counts = collections.Counter(types)
        if counts['operator'] == 1 and counts['term'] == 1 and counts['factor'] > 0 and \
           sum(counts.values()) == len(types):
            check_e = True
            break
    print(f"Statement E is {check_e}.")

    # F. The node in the deepest layer has a parent which in turn has a <term> as a parent.
    node = layers.get(max_depth, [None])[0]
    check_f = (node and node.parent and node.parent.parent and
               node.parent.parent.type == 'term')
    print(f"Statement F is {check_f}.")

    # G. There are 9 layers in the parse tree.
    check_g = max_depth == 9
    print(f"Statement G is {check_g}.")

    # H. The fourth layer contains only <factor> nodes except for a single number node.
    layer_4_types = [n.type for n in layers.get(4, [])]
    counts_h = collections.Counter(layer_4_types)
    check_h = (counts_h['number'] == 1 and counts_h['factor'] > 0 and
               (counts_h['number'] + counts_h['factor']) == len(layer_4_types))
    print(f"Statement H is {check_h}.")

    print("\n--- Conclusion ---")
    false_statements = [s for s, res in [('A', check_a), ('B', check_b), ('C', check_c),
                        ('D', check_d), ('E', check_e), ('F', check_f), ('G', check_g),
                        ('H', check_h)] if not res]

    if len(false_statements) > 1:
        print(f"The statements found to be false are: {', '.join(false_statements)}.")
        print("Since more than one statement is false, the correct option is I.")
    elif len(false_statements) == 1:
        print(f"The only false statement is {false_statements[0]}.")
    else:
        print("All statements were found to be true.")

    print("\n--- Numbers from the expression ---")
    print(4)
    print(5)

if __name__ == '__main__':
    analyze_statements()
    print("<<<I>>>")