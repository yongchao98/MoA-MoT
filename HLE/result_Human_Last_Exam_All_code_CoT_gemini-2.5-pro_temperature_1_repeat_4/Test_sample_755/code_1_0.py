import collections

class Node:
    """A simple class to represent a node in the parse tree."""
    def __init__(self, node_type, value=None):
        self.type = node_type
        self.value = value
        self.children = []
        self.parent = None
        self.depth = 0

    def add_child(self, child_node):
        """Adds a child to this node and sets its parent."""
        self.children.append(child_node)
        child_node.parent = self

    def __repr__(self):
        return f"Node(type='{self.type}', depth={self.depth})"

def build_parse_tree():
    """
    Builds the specific parse tree for the expression 'y + (x + 4) * 5'
    based on the given BNF grammar.
    """
    # Root: <expression> -> <expression> + <term>
    root = Node("<expression>")
    e_y_branch = Node("<expression>")
    op_plus1 = Node("+", "+")
    t_main_right = Node("<term>")
    root.add_child(e_y_branch)
    root.add_child(op_plus1)
    root.add_child(t_main_right)

    # Left branch for 'y': <expression> -> <term> -> <factor> -> name -> y
    t_y = Node("<term>")
    e_y_branch.add_child(t_y)
    f_y = Node("<factor>")
    t_y.add_child(f_y)
    n_y = Node("name")
    f_y.add_child(n_y)
    y_leaf = Node("leaf", "y")
    n_y.add_child(y_leaf)

    # Right branch for '(x + 4) * 5': <term> -> <term> * <factor>
    t_x_plus_4_branch = Node("<term>")
    op_mul = Node("*", "*")
    f_5_branch = Node("<factor>")
    t_main_right.add_child(t_x_plus_4_branch)
    t_main_right.add_child(op_mul)
    t_main_right.add_child(f_5_branch)

    # Sub-branch for '5': <factor> -> number -> 5
    num_5 = Node("number")
    f_5_branch.add_child(num_5)
    five_leaf = Node("leaf", "5")
    num_5.add_child(five_leaf)

    # Sub-branch for '(x + 4)': <term> -> <factor> -> (<expression>)
    f_x_plus_4_paren = Node("<factor>")
    t_x_plus_4_branch.add_child(f_x_plus_4_paren)
    p_left = Node("(", "(")
    e_x_plus_4 = Node("<expression>")
    p_right = Node(")", ")")
    f_x_plus_4_paren.add_child(p_left)
    f_x_plus_4_paren.add_child(e_x_plus_4)
    f_x_plus_4_paren.add_child(p_right)

    # Inside parentheses for 'x + 4': <expression> -> <expression> + <term>
    e_x_branch = Node("<expression>")
    op_plus2 = Node("+", "+")
    t_4_branch = Node("<term>")
    e_x_plus_4.add_child(e_x_branch)
    e_x_plus_4.add_child(op_plus2)
    e_x_plus_4.add_child(t_4_branch)

    # Inner left branch for 'x': <expression> -> <term> -> <factor> -> name -> x
    t_x = Node("<term>")
    e_x_branch.add_child(t_x)
    f_x = Node("<factor>")
    t_x.add_child(f_x)
    n_x = Node("name")
    f_x.add_child(n_x)
    x_leaf = Node("leaf", "x")
    n_x.add_child(x_leaf)

    # Inner right branch for '4': <term> -> <factor> -> number -> 4
    f_4 = Node("<factor>")
    t_4_branch.add_child(f_4)
    num_4 = Node("number")
    f_4.add_child(num_4)
    four_leaf = Node("leaf", "4")
    num_4.add_child(four_leaf)

    return root

def analyze_tree(root):
    """Analyzes the tree to find layers, depth, and all nodes."""
    all_nodes = []
    layers = collections.defaultdict(list)
    q = collections.deque([(root, 1)])
    visited = {root}
    root.depth = 1
    max_depth = 0

    while q:
        curr, depth = q.popleft()
        all_nodes.append(curr)
        layers[depth].append(curr)
        max_depth = max(max_depth, depth)

        for child in curr.children:
            if child not in visited:
                child.depth = depth + 1
                visited.add(child)
                q.append((child, depth + 1))
    
    return all_nodes, layers, max_depth

def evaluate_statements(all_nodes, layers, max_depth):
    """Evaluates statements A-H against the parsed tree."""
    results = {}

    # A. There is at least one <expression> which has a parent that is also an <expression> node.
    check_A = any(n.parent and n.parent.type == "<expression>" for n in all_nodes if n.type == "<expression>")
    results['A'] = check_A

    # B. The deepest number node is in the second to last layer of the tree.
    number_nodes = [n for n in all_nodes if n.type == "number"]
    deepest_num_depth = max(n.depth for n in number_nodes) if number_nodes else -1
    results['B'] = (deepest_num_depth == max_depth - 1)

    # C. There is a name node that appears in a layer which is between ... two layers ... that contain a number node.
    name_depths = {n.depth for n in all_nodes if n.type == "name"}
    num_depths = {n.depth for n in number_nodes}
    check_C = False
    if len(num_depths) >= 2:
        min_d, max_d = min(num_depths), max(num_depths)
        if any(min_d < d < max_d for d in name_depths):
            check_C = True
    results['C'] = check_C

    # D. The deepest layer contains a name with a <factor> as a parent.
    # Interpretation: A leaf node (like 'x') is a "name". Does its parent have type '<factor>'?
    check_D = False
    for node in layers.get(max_depth, []):
        if node.type == 'leaf' and node.value in ['x', 'y']:
             if node.parent and node.parent.type == "<factor>":
                check_D = True
                break
    results['D'] = check_D

    # E. There is a layer that only has <factor> nodes, one operator, and one <term> node.
    check_E = False
    for depth, nodes_in_layer in layers.items():
        types = [n.type for n in nodes_in_layer]
        op_count = sum(1 for t in types if t in ['+', '-', '*', '/'])
        term_count = types.count("<term>")
        factor_count = types.count("<factor>")
        if op_count == 1 and term_count == 1 and factor_count > 0 and (op_count + term_count + factor_count == len(types)):
            check_E = True
            break
    results['E'] = check_E

    # F. The node in the deepest layer has a parent which in turn has a <term> as a parent.
    check_F = False
    for node in layers.get(max_depth, []):
        if node.parent and node.parent.parent and node.parent.parent.type == "<term>":
            check_F = True
            break
    results['F'] = check_F

    # G. There are 9 layers in the parse tree...
    results['G'] = (max_depth == 9)

    # H. The fourth layer contains only <factor> nodes except for a single number node.
    types_in_L4 = [n.type for n in layers.get(4, [])]
    num_count = types_in_L4.count("number")
    factor_count = types_in_L4.count("<factor>")
    results['H'] = (num_count == 1 and (num_count + factor_count == len(types_in_L4)))
    
    return results

def main():
    """Main function to run the analysis and print the result."""
    print("Analyzing the parse tree for: y + (x + 4) * 5\n")
    
    parse_tree_root = build_parse_tree()
    all_nodes, layers, max_depth = analyze_tree(parse_tree_root)
    
    print(f"Parse tree analysis complete. Total layers found: {max_depth}")
    print("-" * 30)

    results = evaluate_statements(all_nodes, layers, max_depth)
    
    false_statements = []
    for statement in sorted(results.keys()):
        is_true = results[statement]
        print(f"Statement {statement}: {is_true}")
        if not is_true:
            false_statements.append(statement)
    
    print("-" * 30)
    print(f"Found {len(false_statements)} false statements: {', '.join(false_statements)}")

    final_answer = ""
    if len(false_statements) > 1:
        final_answer = "I"
        print("\nSince more than one statement is false, the correct option is 'I'.")
    elif len(false_statements) == 1:
        final_answer = false_statements[0]
        print(f"\nSince only statement {final_answer} is false, it is the correct answer.")
    else:
        final_answer = "No false statements found"
        print("\nError: All statements were found to be true.")
    
    print("\nFinal Answer:")
    print(f"<<<{final_answer}>>>")

if __name__ == "__main__":
    main()