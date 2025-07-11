import itertools

def find_identity_and_reason(n, table):
    """Checks for an identity element. Returns (identity, reason_if_failed)."""
    for e in range(n):
        is_identity = True
        for a in range(n):
            if table[e][a] != a or table[a][e] != a:
                is_identity = False
                break
        if is_identity:
            return e, None
    reason = "It fails the Identity Element axiom. No element 'e' exists such that e*a = a and a*e = a for all 'a' in the set."
    return None, reason

def check_inverses_and_reason(n, table, identity):
    """Checks if all elements have an inverse. Returns (bool, reason_if_failed)."""
    if identity is None:
        # This check is technically redundant if called after a failed identity check, but good for completeness.
        return False, "Cannot check for inverses without an identity element."
    for a in range(n):
        has_inverse = False
        for b in range(n):
            if table[a][b] == identity and table[b][a] == identity:
                has_inverse = True
                break
        if not has_inverse:
            reason = f"It fails the Inverse Element axiom. Element '{a}' does not have an inverse with respect to the identity element '{identity}'."
            return False, reason
    return True, None

def check_associativity_and_reason(n, table):
    """Checks for associativity. Returns (bool, reason_if_failed)."""
    for a in range(n):
        for b in range(n):
            for c in range(n):
                # Calculate (a * b) * c
                val1 = table[table[a][b]][c]
                # Calculate a * (b * c)
                val2 = table[a][table[b][c]]
                if val1 != val2:
                    reason = f"It fails the Associativity axiom. For example, with a={a}, b={b}, c={c}: (a*b)*c = {val1} but a*(b*c) = {val2}."
                    return False, reason
    return True, None

def solve():
    """Finds the smallest n for which a non-group structure exists."""
    n = 1
    while True:
        elements = range(n)
        # Generate all possible binary operations for the set {0, ..., n-1}
        op_results_generator = itertools.product(elements, repeat=n * n)

        for op_results in op_results_generator:
            # Reshape the flat list of operation results into a 2D Cayley table
            table = [list(op_results[i * n : (i + 1) * n]) for i in range(n)]

            is_group = True
            failure_reason = ""
            
            # Check axioms in order
            # 1. Identity (Closure is implicit by generation)
            identity, reason = find_identity_and_reason(n, table)
            if identity is None:
                is_group = False
                failure_reason = reason
            else:
                # 2. Inverses
                inverses_ok, reason = check_inverses_and_reason(n, table, identity)
                if not inverses_ok:
                    is_group = False
                    failure_reason = reason
                else:
                    # 3. Associativity
                    associativity_ok, reason = check_associativity_and_reason(n, table)
                    if not associativity_ok:
                        is_group = False
                        failure_reason = reason
            
            # If any axiom failed, we found our answer
            if not is_group:
                print(f"The smallest number n is {n}.")
                print("\nAn example of a set G with n elements and a binary operation `*` that is NOT a group is given by the following operation table:")
                print(f"Let G = {set(elements)}")
                
                header = "  * | " + " ".join(map(str, elements))
                print(header)
                print("----+" + "-" * (len(header) - 3))
                for i, row in enumerate(table):
                    print(f"  {i} | " + " ".join(map(str, row)))
                
                print("\nThis structure is not a group for the following reason:")
                print(failure_reason)

                print("\nThe binary operation is defined by the following equations:")
                for r in range(n):
                    for c in range(n):
                        print(f"{r} * {c} = {table[r][c]}")
                return

        n += 1

solve()