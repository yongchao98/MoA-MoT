def solve_manifold_problem():
    """
    This function programmatically finds the solution to the manifold problem
    by searching for the lexicographically smallest tuple that satisfies the
    conditions derived from the properties of Euler characteristics.
    """
    search_limit = 5  # A small limit is sufficient

    # Genera are non-negative integers, but we must exclude 1.
    possible_genera = [0] + list(range(2, search_limit))
    
    # We determined the minimal length l must be 2.
    # We search for two pairs (a1,b1) and (a2,b2).
    # To find the lexicographically smallest final tuple (a1,b1,a2,b2),
    # we order the search loops and ensure (a1,b1) <= (a2,b2).
    
    for a1 in possible_genera:
        for b1 in possible_genera:
            chi1 = (2 - 2*a1) * (2 - 2*b1)
            
            for a2 in possible_genera:
                # To maintain lexicographical order of pairs, (a1,b1) <= (a2,b2).
                if a2 < a1:
                    continue
                for b2 in possible_genera:
                    if a2 == a1 and b2 < b1:
                        continue

                    chi2 = (2 - 2*a2) * (2 - 2*b2)
                    
                    if chi1 + chi2 == 0:
                        # Since we are iterating in lexicographical order, the first
                        # solution found is the smallest one.
                        
                        # Print the final equation as requested
                        print(f"The determining equation is chi(N) = chi(M({a1},{b1})) + chi(M({a2},{b2})) = 0.")
                        print(f"Numerically: (2 - 2*{a1})*(2 - 2*{b1}) + (2 - 2*{a2})*(2 - 2*{b2}) = {chi1} + {chi2} = 0.")

                        # Print the final tuple in the required format
                        result_tuple = (a1, b1, a2, b2)
                        print(f"The lexicographically least tuple is: ({','.join(map(str, result_tuple))})")
                        return

solve_manifold_problem()