import re

class Node:
    """A simple node for a tree structure."""
    def __init__(self, value, parent=None):
        self.value = value
        self.parent = parent
        self.children = []

    def add_child(self, value):
        """Adds a child node."""
        child_node = Node(value, parent=self)
        self.children.append(child_node)
        return child_node

    def __repr__(self):
        return f"Node({self.value})"

def build_parse_tree():
    """Builds the specific parse tree for y + (x + 4) * 5."""
    root = Node("<expression>") # Layer 1

    # Expansion for y + ...
    expr_left = root.add_child("<expression>")
    root.add_child("+")
    term_right = root.add_child("<term>") # This will be (x+4)*5

    # Left side: y
    term_y = expr_left.add_child("<term>")
    factor_y = term_y.add_child("<factor>")
    factor_y.add_child("name(y)")

    # Right side: (x+4)*5
    term_paren = term_right.add_child("<term>") # This is (x+4)
    term_right.add_child("*")
    factor_5 = term_right.add_child("<factor>")
    factor_5.add_child("number(5)")

    # Expansion for (x+4)
    factor_paren = term_paren.add_child("<factor>")
    factor_paren.add_child("(")
    expr_inner = factor_paren.add_child("<expression>") # This is x+4
    factor_paren.add_child(")")

    # Expansion for x+4
    expr_x = expr_inner.add_child("<expression>")
    expr_inner.add_child("+")
    term_4 = expr_inner.add_child("<term>")

    # Expansion for x
    term_x = expr_x.add_child("<term>")
    factor_x = term_x.add_child("<factor>")
    factor_x.add_child("name(x)")

    # Expansion for 4
    factor_4 = term_4.add_child("<factor>")
    factor_4.add_child("number(4)")

    return root

def get_layers(root_node):
    """Performs a breadth-first search to get nodes layer by layer."""
    if not root_node:
        return []
    layers = []
    queue = [(root_node, 0)]
    visited_nodes = {root_node}

    while queue:
        node, level = queue.pop(0)
        
        # Ensure the layers list is long enough
        while len(layers) <= level:
            layers.append([])
        
        layers[level].append(node)

        for child in node.children:
            if child not in visited_nodes:
                visited_nodes.add(child)
                queue.append((child, level + 1))
    return layers

# --- Functions to check each statement ---

def check_A(layers):
    """A. There is at least one <expression> which has a parent that is also an <expression> node."""
    for layer in layers:
        for node in layer:
            if node.value == "<expression>" and node.parent and node.parent.value == "<expression>":
                return True
    return False

def check_B(layers):
    """B. The deepest number node is in the second to last layer of the tree."""
    deepest_num_layer = -1
    for i, layer in enumerate(layers):
        for node in layer:
            if 'number' in node.value:
                deepest_num_layer = max(deepest_num_layer, i)
    return deepest_num_layer == len(layers) - 2

def check_C(layers):
    """C. There is a name node that appears in a layer which is between two layers that each contain a number node."""
    num_layers = sorted(list({i for i, l in enumerate(layers) for n in l if 'number' in n.value}))
    name_layers = sorted(list({i for i, l in enumerate(layers) for n in l if 'name' in n.value}))
    if len(num_layers) < 2:
        return False
    min_num, max_num = num_layers[0], num_layers[-1]
    for nl in name_layers:
        if min_num < nl < max_num:
            return True
    return False
    
def check_D(layers):
    """D. The deepest layer contains a name with a <factor> as a parent."""
    deepest_layer = layers[-1]
    for node in deepest_layer:
        if 'name' in node.value and node.parent and node.parent.value == "<factor>":
            return True
    return False

def check_E(layers):
    """E. There is a layer that only has <factor> nodes, one operator, and one <term> node."""
    operators = {'+', '-', '*', '/'}
    for layer in layers:
        has_factor = any(node.value == "<factor>" for node in layer)
        op_count = sum(1 for node in layer if node.value in operators)
        term_count = sum(1 for node in layer if node.value == "<term>" for node in layer)
        other_nodes = any(node.value not in ["<factor>", "<term>"] + list(operators) for node in layer)
        
        if has_factor and op_count == 1 and term_count == 1 and not other_nodes:
            return True
    return False

def check_F(layers):
    """F. The node in the deepest layer has a parent which in turn has a <term> as a parent."""
    deepest_layer = layers[-1]
    for node in deepest_layer:
        if node.parent and node.parent.parent and node.parent.value == "<factor>" and node.parent.parent.value == "<term>":
            return True
    return False

def check_G(layers):
    """G. There are 9 layers in the parse tree."""
    return len(layers) == 9

def check_H(layers):
    """H. The fourth layer contains only <factor> nodes except for a single number node."""
    if len(layers) < 4:
        return False
    layer4 = layers[3] # 0-indexed
    num_count = 0
    factor_count = 0
    other_count = 0
    for node in layer4:
        if 'number' in node.value:
            num_count += 1
        elif node.value == '<factor>':
            factor_count += 1
        else:
            other_count += 1
    return num_count == 1 and factor_count > 0 and other_count == 0


def main():
    root = build_parse_tree()
    layers = get_layers(root)

    print("Analysis of statements for the parse tree of 'y + (x + 4) * 5':\n")

    results = {
        'A': check_A(layers),
        'B': check_B(layers),
        'C': check_C(layers),
        'D': check_D(layers),
        'E': check_E(layers),
        'F': check_F(layers),
        'G': check_G(layers),
        'H': check_H(layers),
    }

    false_statement = None
    for statement, is_true in results.items():
        print(f"Statement {statement} is: {is_true}")
        if not is_true:
            false_statement = statement

    print("\n---")
    # Deriving the symbolic equation
    expression = "y + (x + 4) * 5"
    symbolic_expression = re.sub(r'[a-zA-Z_][a-zA-Z0-9_]*', 'name', expression)
    symbolic_expression = re.sub(r'\b\d+\b', 'number', symbolic_expression)
    print(f"The symbolic equation is: {symbolic_expression}")
    print("---\n")
    
    if false_statement:
        print(f"The statement that is NOT true is {false_statement}.")
        print(f"<<<{false_statement}>>>")
    else:
        print("All statements were found to be true.")

if __name__ == "__main__":
    main()
