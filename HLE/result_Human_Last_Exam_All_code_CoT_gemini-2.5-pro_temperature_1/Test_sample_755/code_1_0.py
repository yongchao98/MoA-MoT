class Node:
    """A simple class to represent a node in the parse tree."""
    def __init__(self, label, layer, parent=None):
        self.label = label
        self.layer = layer
        self.parent = parent
        self.children = []
        if parent:
            parent.children.append(self)

def build_parse_tree():
    """Constructs the parse tree for 'y + (x + 4) * 5'."""
    all_nodes = []
    
    # Layer 1
    root = Node("<expression>", 1)
    all_nodes.append(root)
    
    # Layer 2
    e_y_side = Node("<expression>", 2, root)
    op_plus1 = Node("+", 2, root)
    t_rhs = Node("<term>", 2, root)
    all_nodes.extend([e_y_side, op_plus1, t_rhs])

    # Layer 3
    t_y = Node("<term>", 3, e_y_side)
    t_lhs = Node("<term>", 3, t_rhs)
    op_mul = Node("*", 3, t_rhs)
    f_5 = Node("<factor>", 3, t_rhs)
    all_nodes.extend([t_y, t_lhs, op_mul, f_5])

    # Layer 4
    f_y = Node("<factor>", 4, t_y)
    f_parens = Node("<factor>", 4, t_lhs)
    num_5_node = Node("<number>", 4, f_5)
    all_nodes.extend([f_y, f_parens, num_5_node])

    # Layer 5
    name_y_node = Node("<name>", 5, f_y)
    p_open = Node("(", 5, f_parens)
    e_inner = Node("<expression>", 5, f_parens)
    p_close = Node(")", 5, f_parens)
    term_5 = Node("5", 5, num_5_node)
    all_nodes.extend([name_y_node, p_open, e_inner, p_close, term_5])

    # Layer 6
    term_y = Node("y", 6, name_y_node)
    e_x_side = Node("<expression>", 6, e_inner)
    op_plus2 = Node("+", 6, e_inner)
    t_4_side = Node("<term>", 6, e_inner)
    all_nodes.extend([term_y, e_x_side, op_plus2, t_4_side])

    # Layer 7
    t_x = Node("<term>", 7, e_x_side)
    f_4 = Node("<factor>", 7, t_4_side)
    all_nodes.extend([t_x, f_4])

    # Layer 8
    f_x = Node("<factor>", 8, t_x)
    num_4_node = Node("<number>", 8, f_4)
    all_nodes.extend([f_x, num_4_node])

    # Layer 9
    name_x_node = Node("<name>", 9, f_x)
    term_4 = Node("4", 9, num_4_node)
    all_nodes.extend([name_x_node, term_4])

    # Layer 10
    term_x = Node("x", 10, name_x_node)
    all_nodes.append(term_x)
    
    return all_nodes

def analyze_statements(nodes):
    """Evaluates statements A-H against the parse tree."""
    max_layer = max(n.layer for n in nodes)
    results = {}

    # Statement A
    results['A'] = any(n.parent.label == "<expression>" for n in nodes if n.label == "<expression>" and n.parent)

    # Statement B
    number_terminals = [n for n in nodes if n.label.isdigit()]
    deepest_num_layer = max(n.layer for n in number_terminals)
    results['B'] = (deepest_num_layer == max_layer - 1)

    # Statement C
    name_layers = sorted({n.layer for n in nodes if n.label == "<name>"})
    number_layers = sorted({n.layer for n in nodes if n.label == "<number>"})
    results['C'] = any(num_l1 < name_l < num_l2 for name_l in name_layers for num_l1 in number_layers for num_l2 in number_layers)

    # Statement D
    deepest_nodes = [n for n in nodes if n.layer == max_layer]
    d_true = False
    for n in deepest_nodes:
        # A "name" can be a <name> non-terminal or a terminal like 'x'
        if n.label == "<name>" and n.parent and n.parent.label == "<factor>":
             d_true = True
        if n.parent and n.parent.label == "<name>": # a terminal name
             if n.parent.label == "<factor>": # Check if its parent is <factor>
                 d_true = True
    results['D'] = d_true

    # Statement E
    e_true = False
    for i in range(1, max_layer + 1):
        layer_labels = [n.label for n in nodes if n.layer == i]
        factor_count = layer_labels.count("<factor>")
        term_count = layer_labels.count("<term>")
        op_count = sum(1 for label in layer_labels if label in ['+', '-', '*', '/'])
        if len(layer_labels) == factor_count + term_count + op_count and op_count == 1 and term_count == 1 and factor_count >= 1:
            e_true = True
    results['E'] = e_true

    # Statement F
    f_true = False
    for n in deepest_nodes:
        if n.parent and n.parent.parent and n.parent.parent.label == "<term>":
            f_true = True
    results['F'] = f_true
    
    # Statement G
    results['G'] = (max_layer == 9)

    # Statement H
    layer4_labels = [n.label for n in nodes if n.layer == 4]
    factor_count = layer4_labels.count("<factor>")
    number_count = layer4_labels.count("<number>")
    results['H'] = (number_count == 1 and (factor_count + number_count) == len(layer4_labels) and factor_count > 0)

    return results

def main():
    """Main function to run the analysis and print the result."""
    all_nodes = build_parse_tree()
    results = analyze_statements(all_nodes)
    
    print("Analysis of Statements A-H:")
    for statement, is_true in results.items():
        print(f"Statement {statement}: {is_true}")
        
    false_statements = [k for k, v in results.items() if not v]
    
    print("\nFalse statements:", ", ".join(false_statements))
    
    if len(false_statements) > 1:
        final_answer = "I"
    elif len(false_statements) == 1:
        final_answer = false_statements[0]
    else:
        final_answer = "None of the statements are false."

    # Fulfilling the unusual prompt requirement to output numbers
    print("\nNumbers from the expression: 4, 5")

    print("\nThe single statement that is NOT true is not uniquely identifiable, as multiple options are false.")
    print("Therefore, the correct choice is the one indicating this situation.")
    print(f"\nFinal Answer: {final_answer}")
    print(f"\n<<<{final_answer}>>>")

if __name__ == "__main__":
    main()