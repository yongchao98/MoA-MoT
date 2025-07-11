import math

def lcm(a, b):
    """Computes the least common multiple of two integers."""
    return abs(a*b) // math.gcd(a, b) if a != 0 and b != 0 else 0

def group_exponent(orders):
    """Computes the exponent of a group given the orders of its elements."""
    if not orders:
        return 1
    result = orders[0]
    for i in range(1, len(orders)):
        result = lcm(result, orders[i])
    return result

def solve_graph_problem():
    """
    Solves the problem by identifying groups, finding exponents, and summing them by column.
    """
    # Step 1 & 2: Define groups and their exponents
    # Exponent is the lcm of the orders of all elements in the group.
    exponents = {
        'PSL(2,4)': 30,  # isomorphic to A5, orders are 1, 2, 3, 5
        'Z4xZ4': 4,      # max order of element is 4
        'D8': 8,         # dihedral of order 16, has an element of order 8
        'S4': 12,        # orders are 1, 2, 3, 4
        'A4': 6,         # orders are 1, 2, 3
        'D3': 6,         # orders are 1, 2, 3
        'Z3xZ3': 3,      # max order of element is 3
        'Z2_3': 2        # all non-identity elements have order 2
    }

    # Step 3: Map visualizations in each column to their identified group
    # Based on vertex counting and resolving inconsistencies as explained in the text.
    # V6 is assumed to represent Z3xZ3 with identity omitted.
    # V11 is assumed to represent Z4xZ4 despite its ambiguous appearance.
    column_groups = {
        1: ['PSL(2,4)', 'Z4xZ4', 'PSL(2,4)', 'Z2_3'],      # V1, V5, V9, V13
        2: ['D8', 'Z3xZ3', 'Z3xZ3', 'S4'],                # V2, V6, V10, V14
        3: ['A4', 'D3', 'Z4xZ4', 'D8'],                   # V3, V7, V11, V15
        4: ['S4', 'Z2_3', 'D3', 'A4']                     # V4, V8, V12, V16
    }
    
    # Step 4: Calculate and print the sum of exponents for each column
    column_sums = []
    for i in range(1, 5):
        groups_in_col = column_groups[i]
        exp_values = [exponents[group] for group in groups_in_col]
        total_sum = sum(exp_values)
        column_sums.append(total_sum)
        
        # Print the equation for each sum
        exp_str = ' + '.join(map(str, exp_values))
        print(f"S{i} = {exp_str} = {total_sum}")
        
    print("\nThe ordered list of the four column sums is:")
    print(column_sums)

solve_graph_problem()