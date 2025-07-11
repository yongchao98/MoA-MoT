import sys

def solve():
    """
    Solves the problem by searching for the lexicographically smallest tuple.
    """
    
    # The number of manifolds, l, must be 3.
    l = 3
    target_chi_sum = 2 * (l - 1)

    # Search limit for genus a and b. The solution is expected to use small numbers.
    search_limit = 5

    # Generate a list of possible pairs (a,b) and their Euler characteristics chi.
    # The pairs are generated in lexicographical order.
    # a,b are non-negative integers, but not equal to 1.
    pairs = []
    for a in range(search_limit):
        if a == 1:
            continue
        for b in range(search_limit):
            if b == 1:
                continue
            # Euler characteristic of M(a) x M(b) is 4*(1-a)*(1-b)
            chi = 4 * (1 - a) * (1 - b)
            pairs.append(((a, b), chi))

    # Find the first combination of 3 pairs (p1, p2, p3) that satisfies the sum.
    # The loops are structured to find the lexicographically smallest tuple
    # (a1, b1, a2, b2, a3, b3) by ensuring p1 <= p2 <= p3.
    for i1 in range(len(pairs)):
        p1, chi1 = pairs[i1]
        for i2 in range(i1, len(pairs)):
            p2, chi2 = pairs[i2]
            for i3 in range(i2, len(pairs)):
                p3, chi3 = pairs[i3]

                if chi1 + chi2 + chi3 == target_chi_sum:
                    # Found the lexicographically smallest solution.
                    result_tuple = p1 + p2 + p3
                    
                    # Print the final equation with the numbers.
                    # M_sum = M(a1,b1) # M(a2,b2) # M(a3,b3)
                    # chi(M_sum) = chi1 + chi2 + chi3 - 2*(l-1) = 0
                    # chi1 + chi2 + chi3 = 4
                    print(f"The Euler characteristics of the three manifolds are: {chi1}, {chi2}, {chi3}")
                    print(f"The equation for the Euler characteristic of the connect-sum is:")
                    print(f"{chi1} + {chi2} + {chi3} = {target_chi_sum}")

                    # Print the final tuple as requested.
                    flat_tuple_str = "(" + ",".join(map(str, result_tuple)) + ")"
                    print(f"\nThe lexicographically least tuple is:")
                    print(flat_tuple_str)
                    return

if __name__ == '__main__':
    solve()
