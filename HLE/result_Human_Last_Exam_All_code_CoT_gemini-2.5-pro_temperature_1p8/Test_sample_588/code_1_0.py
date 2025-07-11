import sys

def solve():
    """
    Solves the topology problem to find the lexicographically least tuple
    based on the analysis of Euler characteristics.
    """

    def chi(a, b):
        """Calculates the Euler characteristic of M(a,b)."""
        return (2 - 2 * a) * (2 - 2 * b)

    # A search limit sufficient for finding the smallest pairs
    max_search_val = 5
    pairs = []
    # Generate candidate (a,b) pairs (a!=1, b!=1), sorted lexicographically.
    # We assume a <= b since M(a,b) is symmetric w.r.t a and b.
    for a in range(max_search_val):
        if a == 1:
            continue
        for b in range(a, max_search_val):
            if b == 1:
                continue
            pairs.append((a, b))

    # From our analysis, the minimal number of manifolds is l=3.
    # The target sum for the Euler characteristics is 2*(3-1) = 4.
    l = 3
    target_sum = 4

    # Search for the lexicographically smallest set of three pairs {p1, p2, p3}
    # whose chi values sum to the target. The loops are nested to find the
    # first possible combination in lexicographical order.
    for i in range(len(pairs)):
        p1 = pairs[i]
        chi1 = chi(p1[0], p1[1])
        # The loops for j and k start from i and j respectively to ensure the
        # set of pairs (p1, p2, p3) is sorted, guaranteeing the final tuple is minimal.
        for j in range(i, len(pairs)):
            p2 = pairs[j]
            chi2 = chi(p2[0], p2[1])
            for k in range(j, len(pairs)):
                p3 = pairs[k]
                chi3 = chi(p3[0], p3[1])

                if chi1 + chi2 + chi3 == target_sum:
                    # The first solution found is guaranteed to be the minimal one.
                    result_pairs = [p1, p2, p3]
                    chi_values = [chi1, chi2, chi3]
                    
                    print(f"Minimal number of manifolds, l, is {l}.")
                    print(f"The governing equation is Sum(χ_i) = 2*(l-1), which for l=3 is {target_sum}.")
                    print("The lexicographically smallest set of non-full manifolds M(a,b) that solve the equation are:")
                    print(f"  M({result_pairs[0][0]},{result_pairs[0][1]}) with χ = {chi_values[0]}")
                    print(f"  M({result_pairs[1][0]},{result_pairs[1][1]}) with χ = {chi_values[1]}")
                    print(f"  M({result_pairs[2][0]},{result_pairs[2][1]}) with χ = {chi_values[2]}")
                    
                    print("\nThe equation with these manifolds is:")
                    print(f"χ(M({p1[0]},{p1[1]})) + χ(M({p2[0]},{p2[1]})) + χ(M({p3[0]},{p3[1]})) = {chi1} + {chi2} + {chi3} = {target_sum}")
                    
                    final_tuple = (p1[0], p1[1], p2[0], p2[1], p3[0], p3[1])
                    print("\nThe final lexicographically least tuple (a1,b1,a2,b2,a3,b3) is:")
                    print(f"({','.join(map(str, final_tuple))})")

                    return final_tuple

# Execute the solver and format the final answer.
final_tuple = solve()
flat_tuple_str = f"({','.join(map(str, final_tuple))})"
print(f"<<<{flat_tuple_str}>>>")
