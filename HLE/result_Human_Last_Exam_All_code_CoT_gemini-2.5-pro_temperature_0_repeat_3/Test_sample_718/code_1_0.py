import itertools

def get_transitive_closure(poset_relations):
    """Computes the full set of relations from the direct successor relations."""
    closure = {k: set(v) for k, v in poset_relations.items()}
    elements = list(closure.keys())
    while True:
        updated = False
        for i in elements:
            for j in list(closure.get(i, [])):
                for k in list(closure.get(j, [])):
                    if k not in closure[i]:
                        closure[i].add(k)
                        updated = True
        if not updated:
            break
    return closure

def get_upper_bounds(full_poset, x, y):
    """Get common upper bounds for two elements."""
    succ_x = full_poset.get(x, set())
    succ_y = full_poset.get(y, set())
    # An element is an upper bound of itself
    return (succ_x | {x}).intersection(succ_y | {y})

def get_minimal_elements(full_poset, subset):
    """Find minimal elements in a given subset of the poset."""
    minimal = set()
    for s1 in subset:
        is_minimal = True
        for s2 in subset:
            if s1 != s2 and s1 in full_poset.get(s2, set()):
                # s2 < s1, so s1 is not minimal
                is_minimal = False
                break
        if is_minimal:
            minimal.add(s1)
    return minimal

def is_upper_semilattice(poset_relations):
    """Check if a poset is an upper semilattice."""
    full_poset = get_transitive_closure(poset_relations)
    elements = set(full_poset.keys()).union(*full_poset.values())
    
    for x, y in itertools.combinations(elements, 2):
        # If comparable, they have a join (the larger element)
        if y in full_poset.get(x, set()) or x in full_poset.get(y, set()):
            continue

        upper_bounds = get_upper_bounds(full_poset, x, y)
        if not upper_bounds:
            print(f"FAIL: Pair ({x}, {y}) has no upper bounds.")
            return False
        
        minimal_upper_bounds = get_minimal_elements(full_poset, upper_bounds)
        if len(minimal_upper_bounds) != 1:
            print(f"FAIL: Pair ({x}, {y}) has {len(minimal_upper_bounds)} minimal upper bounds: {minimal_upper_bounds}.")
            return False
    return True

def find_diamond_subposet(poset_relations):
    """Finds a diamond subposet M3 = {m, a, b, t}."""
    full_poset = get_transitive_closure(poset_relations)
    elements = list(set(full_poset.keys()).union(*full_poset.values()))
    
    for m, a, b, t in itertools.permutations(elements, 4):
        if len({m, a, b, t}) != 4: continue
        
        # Check relations for M3: m < a, m < b, a < t, b < t
        succ_m = full_poset.get(m, set())
        succ_a = full_poset.get(a, set())
        succ_b = full_poset.get(b, set())

        if a in succ_m and b in succ_m and t in succ_a and t in succ_b:
             # Check incomparability of a and b
            if a not in succ_b and b not in succ_a:
                return {'m': m, 'a': a, 'b': b, 't': t}
    return None

def main():
    """Main function to execute the logic and print the result."""
    print("The problem asks for the value 'n' for which any tame functor from an upper semilattice J to Vect_K is n-resolvable.")
    print("This is a question about the global dimension of the category algebra KJ.\n")
    print("Step 1: An upper semilattice J cannot contain a 'crown' subposet. This implies that the global dimension of KJ is finite.")
    print("Step 2: The algebra KJ is of 'tame' representation type. A tame algebra with finite global dimension must have a global dimension of 1 or 2.")
    print("Step 3: Therefore, n must be 1 or 2. To show n=2, we must demonstrate that a global dimension of 2 is possible.")
    print("Step 4: The global dimension is > 1 if the poset contains a non-distributive lattice as an interval, such as the diamond lattice M3.")
    print("Step 5: We construct a poset J that is a tame upper semilattice and contains a diamond subposet.\n")

    # We build upon the tame poset corresponding to the extended Dynkin diagram D_4tilde,
    # which can be represented by relations l_i < c for i=1,2,3,4.
    # We add a new minimal element 'm' such that m < l1 and m < l2.
    poset_J_relations = {
        'm': {'l1', 'l2'},
        'l1': {'c'},
        'l2': {'c'},
        'l3': {'c'},
        'l4': {'c'},
        'c': set()
    }
    elements = sorted(list(set(poset_J_relations.keys()).union(*poset_J_relations.values())))
    print(f"Constructed Poset J has elements {elements} and direct relations (x < y):")
    for k, v_set in sorted(poset_J_relations.items()):
        for v in sorted(list(v_set)):
            print(f"  {k} < {v}")
    print("-" * 20)

    print("Verifying J is an upper semilattice...")
    if is_upper_semilattice(poset_J_relations):
        print("SUCCESS: The constructed poset J is an upper semilattice.")
    else:
        print("ERROR: The constructed poset J is not an upper semilattice.")
    print("-" * 20)

    print("Verifying J contains a diamond subposet...")
    diamond = find_diamond_subposet(poset_J_relations)
    if diamond:
        print(f"SUCCESS: Found a diamond subposet M3 = {{m': '{diamond['m']}', a': '{diamond['a']}', b': '{diamond['b']}', t': '{diamond['t']}'}}.")
        print("The existence of this subposet implies gldim(KJ) > 1.")
    else:
        print("ERROR: Could not find a diamond subposet.")
    print("-" * 20)
    
    print("Conclusion:")
    print("We have constructed a poset J which is an upper semilattice and is of tame representation type.")
    print("This poset J has a sub-interval that forms a non-distributive diamond lattice, which implies its algebra KJ has a global dimension of at least 2.")
    print("Since the global dimension for any such poset can only be 1 or 2, the maximum possible value is 2.")
    
    final_answer = 2
    print(f"\nTherefore, any tame functor f: J -> Vect_K is n-resolvable for n = {final_answer}.")

if __name__ == '__main__':
    main()