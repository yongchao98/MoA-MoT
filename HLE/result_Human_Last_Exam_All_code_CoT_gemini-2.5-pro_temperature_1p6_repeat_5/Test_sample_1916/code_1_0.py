import itertools

def check_associativity(table):
    """Checks if the multiplication defined by the table is associative."""
    n = len(table)
    for x in range(n):
        for y in range(n):
            for z in range(n):
                # Check if (x*y)*z == x*(y*z)
                left_side = table[table[x][y]][z]
                right_side = table[x][table[y][z]]
                if left_side != right_side:
                    return False
    return True

def are_isomorphic(table1, table2):
    """
    Checks if two monoids are isomorphic.
    An isomorphism must map identity to identity (0 -> 0).
    We only need to check the two permutations of the non-identity elements {1, 2}.
    """
    # Permutation 1: f(1)=1, f(2)=2
    f1 = [0, 1, 2]
    is_iso1 = True
    for i in range(3):
        for j in range(3):
            # Check if f(i*j) == f(i)*f(j)
            if f1[table1[i][j]] != table2[f1[i]][f1[j]]:
                is_iso1 = False
                break
        if not is_iso1:
            break
    if is_iso1:
        return True

    # Permutation 2: f(1)=2, f(2)=1
    f2 = [0, 2, 1]
    is_iso2 = True
    for i in range(3):
        for j in range(3):
            # Check if f(i*j) == f(i)*f(j)
            if f2[table1[i][j]] != table2[f2[i]][f2[j]]:
                is_iso2 = False
                break
        if not is_iso2:
            break
    if is_iso2:
        return True

    return False

def main():
    """
    Main function to find the number of non-isomorphic monoids of order 3.
    """
    elements = [0, 1, 2]  # Representing {e, a, b}
    
    # The multiplication table is determined by the 4 products of non-identity elements
    num_potential_structures = len(elements)**4
    print(f"Total potential structures to check: 3^4 = {num_potential_structures}")

    valid_monoids = []
    # Iterate through all 81 possibilities for the 2x2 sub-table
    for aa, ab, ba, bb in itertools.product(elements, repeat=4):
        # Construct the full 3x3 Cayley table
        table = [
            [0, 1, 2],    # Row for identity element 'e' (0)
            [1, aa, ab],  # Row for element 'a' (1)
            [2, ba, bb]   # Row for element 'b' (2)
        ]
        
        if check_associativity(table):
            valid_monoids.append(table)
            
    print(f"Number of associative structures found: {len(valid_monoids)}")

    # Filter the list of valid monoids to find unique isomorphism classes
    unique_monoids = []
    for monoid_table in valid_monoids:
        is_new = True
        for unique_table in unique_monoids:
            if are_isomorphic(monoid_table, unique_table):
                is_new = False
                break
        if is_new:
            unique_monoids.append(monoid_table)
            
    print(f"Number of non-isomorphic classes found: {len(unique_monoids)}")

if __name__ == '__main__':
    main()