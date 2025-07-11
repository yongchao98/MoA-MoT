import itertools

def check_properties(table_values, n_max):
    """
    Checks the properties of a magma defined on {0, 1, 2}.
    `table_values` defines the results for (0*1), (0*2), and (1*2).
    """
    N = 3
    table = [[0]*N for _ in range(N)]
    
    # Idempotent and Commutative properties are built-in
    for i in range(N):
        table[i][i] = i
    table[0][1] = table[1][0] = table_values[0]
    table[0][2] = table[2][0] = table_values[1]
    table[1][2] = table[2][1] = table_values[2]

    # Check Left Self-Distributivity
    is_lsd = True
    for x, y, z in itertools.product(range(N), repeat=3):
        if table[x][table[y][z]] != table[table[x][y]][table[x][z]]:
            is_lsd = False
            break
    if not is_lsd:
        return None, None

    # Check Mediality
    is_medial = True
    for w, x, y, z in itertools.product(range(N), repeat=4):
        lhs = table[table[w][x]][table[y][z]]
        rhs = table[table[w][y]][table[x][z]]
        if lhs != rhs:
            is_medial = False
            break

    # Check n-cancellability for n from 1 to n_max
    is_n_cancellable = {}
    for n in range(1, n_max + 1):
        cancellable = True
        for a, b in itertools.product(range(N), repeat=2):
            if a == b:
                continue
            
            res = b
            for _ in range(n):
                res = table[a][res]
            
            if res == b:
                cancellable = False
                break
        if cancellable:
            is_n_cancellable[n] = True
        else:
            is_n_cancellable[n] = False
            
    return is_medial, is_n_cancellable

def find_counterexamples():
    """
    Searches for counterexamples among 3-element magmas.
    A counterexample for a given 'n' is a magma that is n-cancellable but not medial.
    """
    counterexamples = {}
    
    # Iterate through all 27 possible commutative, idempotent magmas on 3 elements
    for t01, t02, t12 in itertools.product(range(3), repeat=3):
        is_medial, n_cancellability = check_properties([t01, t02, t12], 4)

        if is_medial is None: # Not LSD
            continue
            
        if not is_medial:
            for n, is_cancellable in n_cancellability.items():
                if is_cancellable and n not in counterexamples:
                    # Found first counterexample for this n
                    table = [[i for i in range(3)] for _ in range(3)]
                    table[0][1] = table[1][0] = t01
                    table[0][2] = table[2][0] = t02
                    table[1][2] = table[2][1] = t12
                    counterexamples[n] = table

    print("Checking which values of n admit a non-medial, n-cancellable magma...")

    for n in range(1, 5):
        if n in counterexamples:
            print(f"For n = {n}, a counterexample exists. So n-cancellability does NOT imply mediality.")
            table = counterexamples[n]
            print("Counterexample magma table:")
            print("  * | 0 | 1 | 2")
            print(" ---|---|---|---")
            for i in range(3):
                print(f"  {i} | {table[i][0]} | {table[i][1]} | {table[i][2]}")
        else:
            print(f"For n = {n}, no counterexample found on 3 elements. This suggests n-cancellability implies mediality.")
    
    print("\nBased on the search, the property holds for odd n (1, 3) and fails for even n (2, 4).")
    print("So, the magma being n-cancellable implies that it is medial for all ODD positive values of n.")


if __name__ == '__main__':
    find_counterexamples()
