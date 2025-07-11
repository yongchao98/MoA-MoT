def check_structure_properties(n, elements, table):
    """
    Checks if a given algebraic structure is a group and prints the results.
    - n: The number of elements.
    - elements: A list of the elements (e.g., [0, 1]).
    - table: The Cayley table (as a list of lists) for the operation.
    """
    print(f"--- Checking for a structure with n = {n} ---")
    print("Elements: ", elements)
    print("Cayley Table:")
    header = "  | " + " ".join(map(str, elements))
    print(header)
    print("-" * len(header))
    for i, row in enumerate(table):
        print(f"{elements[i]} | " + " ".join(map(str, row)))
    
    # 1. Check Associativity
    is_associative = True
    for i in range(n):
        for j in range(n):
            for k in range(n):
                # (i.j).k vs i.(j.k)
                # Note: table indices correspond to element indices
                lhs = table[table[i][j]][k]
                rhs = table[i][table[j][k]]
                if lhs != rhs:
                    is_associative = False
                    break
            if not is_associative: break
        if not is_associative: break
    print(f"\n1. Is associative: {is_associative}")

    # 2. Check for an Identity Element
    identity_element = None
    for e_idx, e in enumerate(elements):
        is_identity = True
        for a_idx, a in enumerate(elements):
            if table[e_idx][a_idx] != a or table[a_idx][e_idx] != a:
                is_identity = False
                break
        if is_identity:
            identity_element = e
            break
    
    print(f"2. Has an identity element: {identity_element is not None}")
    if identity_element is not None:
        print(f"   - Identity element is: {identity_element}")

    # 3. Check for Inverse Elements
    has_inverses = True
    if identity_element is None:
        has_inverses = False
        print("3. Has an inverse for each element: False (since no identity element exists)")
    else:
        id_idx = elements.index(identity_element)
        for i, elem_a in enumerate(elements):
            found_inverse = False
            for j, elem_b in enumerate(elements):
                if table[i][j] == identity_element and table[j][i] == identity_element:
                    found_inverse = True
                    break
            if not found_inverse:
                has_inverses = False
                break
        print(f"3. Has an inverse for each element: {has_inverses}")

    is_group = is_associative and (identity_element is not None) and has_inverses
    print(f"\nConclusion: The structure IS a group. -> {is_group}")
    print("-" * (len(header) + 2))
    return is_group

# Case n=1
n1_elements = [0]
# The only possible operation: 0.0 = 0
n1_table = [[0]]
check_structure_properties(1, n1_elements, n1_table)

# Case n=2
n2_elements = [0, 1]
# An example operation that is NOT a group: a.b = a
n2_table = [[0, 0], [1, 1]]
is_n2_group = check_structure_properties(2, n2_elements, n2_table)

if not is_n2_group:
    print("\nSince n=1 must form a group, and we found a non-group for n=2,")
    print("the smallest number n is 2.")
    final_answer = 2
    print(f"The final answer is {final_answer}")

<<<2>>>