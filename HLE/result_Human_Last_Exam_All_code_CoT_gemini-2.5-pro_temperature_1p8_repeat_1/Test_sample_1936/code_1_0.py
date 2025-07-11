import itertools

def check_group_axioms(n, table):
    """
    Checks if a given operation table defines a group on the set {0, 1, ..., n-1}.

    Args:
        n: The number of elements in the set.
        table: An n x n list of lists representing the operation table.

    Returns:
        A tuple (is_group, reason), where is_group is a boolean and
        reason is a string explaining why it is or is not a group.
    """
    # Axiom 1: Closure
    # Closure is guaranteed by the way the tables are generated, as all
    # results are in {0, 1, ..., n-1}.

    # Axiom 3: Identity Element
    identity_element = -1
    for e in range(n):
        is_identity_candidate = True
        for i in range(n):
            if table[e][i] != i or table[i][e] != i:
                is_identity_candidate = False
                break
        if is_identity_candidate:
            identity_element = e
            break

    if identity_element == -1:
        return False, "it lacks an identity element"

    # Axiom 4: Inverse Element
    for i in range(n):
        has_inverse = False
        for j in range(n):
            if table[i][j] == identity_element and table[j][i] == identity_element:
                has_inverse = True
                break
        if not has_inverse:
            return False, f"element {i} does not have an inverse"

    # Axiom 2: Associativity
    for i in range(n):
        for j in range(n):
            for k in range(n):
                # (i . j) . k
                lhs = table[table[i][j]][k]
                # i . (j . k)
                rhs = table[i][table[j][k]]
                if lhs != rhs:
                    # Found a failure, provide the specific numbers in the "equation"
                    reason = f"it is not associative. For example: ({i} . {j}) . {k} = {lhs}, but {i} . ({j} . {k}) = {rhs}"
                    return False, reason

    return True, "it satisfies all group axioms"


def find_smallest_n_for_non_group():
    """
    Finds the smallest n for which a non-group structure can be defined.
    """
    n = 1
    while True:
        print(f"--- Checking for n = {n} ---")
        elements = range(n)
        
        # An iterator for all possible values for the n*n cells of the table
        # There are n^(n*n) possible binary operations.
        possible_tables = itertools.product(elements, repeat=n*n)

        found_non_group = False
        for table_flat in possible_tables:
            # Reshape the flat tuple into an n x n table
            table = [list(table_flat[i*n : (i+1)*n]) for i in range(n)]

            is_group, reason = check_group_axioms(n, table)
            
            if not is_group:
                print(f"Found a structure with {n} element(s) that is NOT a group.")
                print(f"The reason is that {reason}.")
                print("\nSet G = {" + ", ".join(map(str, elements)) + "}")
                print("Example of a failing binary operation table:")
                header = " . | " + " ".join(map(str, elements))
                print(header)
                print("-" * len(header))
                for i in range(n):
                    row_str = f" {i} | " + " ".join(map(str, table[i]))
                    print(row_str)
                
                print(f"\nThus, the smallest number n is {n}.")
                found_non_group = True
                break
        
        if found_non_group:
            return n
            
        print(f"All {n**(n*n)} possible binary operations for n={n} form a group.")
        n += 1

if __name__ == '__main__':
    find_smallest_n_for_non_group()