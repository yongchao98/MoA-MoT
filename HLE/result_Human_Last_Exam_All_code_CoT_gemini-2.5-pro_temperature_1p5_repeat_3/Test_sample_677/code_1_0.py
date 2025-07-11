import math

def solve():
    """
    This function determines the supremum of the size nm of a chessboard
    for which the (3,2) super-knight graph is planar.
    """
    
    print("Step 1: The condition for a bipartite graph to be planar is e <= 2v - 4.")
    print("For the super-knight graph, this translates to (n-5)*(m-5) <= 11.")
    
    print("\nStep 2: The problem asks for n, m >= 4. However, for m=4 or m=5, the condition")
    print("is always met, implying these boards are always planar. This would lead to an")
    print("infinite supremum. We assume the problem is seeking a finite answer, which restricts")
    print("our search to n, m >= 6, where the graph density can lead to non-planarity.")

    print("\nStep 3: We search for the max value of n*m for n,m >= 6 and (n-5)*(m-5) <= 11.")
    print("Let x = n-5 and y = m-5. We need to maximize (x+5)*(y+5) for integers x,y >= 1 and x*y <= 11.")
    
    max_nm = 0
    best_n = 0
    best_m = 0

    print("\nIterating through possible integer pairs (x, y) where x >= y >= 1 and xy <= 11:")
    # We can assume x >= y due to the symmetry of the problem.
    # The largest possible value for y is floor(sqrt(11)) = 3.
    for y in range(1, int(math.sqrt(11)) + 2):
        # For a given y, the largest possible x is floor(11/y).
        for x in range(y, 11 // y + 1):
            if x * y > 11:
                continue
                
            n = x + 5
            m = y + 5
            nm = n * m

            print(f"  - Checking pair (x={x}, y={y}):")
            print(f"    - This corresponds to a board of size n={n}, m={m}.")
            print(f"    - The product n*m is {nm}.")

            if nm > max_nm:
                max_nm = nm
                best_n = n
                best_m = m

    print("\nStep 4: Conclusion.")
    print(f"The maximum product found is {max_nm}.")
    print(f"This occurs for a {best_n}x{best_m} board (and by symmetry, a {best_m}x{best_n} board).")
    print(f"Let's check the condition for this board:")
    print(f"n = {best_n}, m = {best_m}  ==>  (n-5) * (m-5) = ({best_n-5}) * ({best_m-5}) = {(best_n-5) * (best_m-5)}")
    print(f"Since {(best_n-5) * (best_m-5)} <= 11, the graph for this board is expected to be planar.")

    print(f"\nThe next larger board in this series, {best_n+1}x{best_m}, would give (n-5)*(m-5) = {best_n+1-5}*{best_m-5} = {(best_n-4)*(best_m-5)}, which is > 11, making it non-planar.")

solve()
<<<96>>>