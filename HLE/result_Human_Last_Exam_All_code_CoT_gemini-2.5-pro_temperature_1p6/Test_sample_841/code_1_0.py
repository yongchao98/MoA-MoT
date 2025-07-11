import sys

def solve_unique_length_factorizations():
    """
    Solves the problem by counting the number of specified rings
    that are Half-Factorial Domains (HFDs).
    """

    # Step 1: Rings that are integrally closed (Maximal Orders O_K)
    # These are HFDs if their class number h(-d) is 1 or 2.
    
    # List of square-free d > 0 where h(-d) = 1 (from Stark-Heegner theorem)
    d_h1 = [1, 2, 3, 7, 11, 19, 43, 67, 163]
    
    # List of square-free d > 0 where h(-d) = 2 (from Baker, Stark)
    d_h2 = [5, 6, 10, 13, 15, 22, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427]

    print("--- Part 1: Integrally Closed Rings (Maximal Orders) ---")
    print("A ring of integers O(Q(sqrt(-d))) has factorizations of unique length if its class number h(-d) is 1 or 2.")
    
    count_h1 = len(d_h1)
    print(f"Number of rings with h(-d)=1 (UFDs): {count_h1}")
    print(f"Values of d: {d_h1}")
    
    count_h2 = len(d_h2)
    print(f"Number of rings with h(-d)=2: {count_h2}")
    print(f"Values of d: {d_h2}")

    count_maximal_HFDs = count_h1 + count_h2
    print(f"\nTotal count for integrally closed rings = {count_h1} + {count_h2} = {count_maximal_HFDs}")
    
    # Step 2: Rings that are not integrally closed (Non-Maximal Orders)
    # These are rings R = Z[sqrt(-d)] where d is square-free and d = 3 (mod 4).
    # These rings are HFDs if their order class number h(R) is 1 or 2.
    # The formula for the order class number h(R) simplifies based on d mod 8:
    # - If d = 3: h(Z[sqrt(-3)]) = 1. This is an HFD.
    # - If d > 3 and d = 3 (mod 8): h(R) = h(-d). R is an HFD if h(-d) <= 2.
    # - If d = 7 (mod 8): h(R) = 3*h(-d). R is never an HFD as 3*h(-d) > 2.

    print("\n--- Part 2: Rings Not Integrally Closed ---")
    print("These are rings Z[sqrt(-d)] for square-free d = 3 (mod 4).")
    print("They have unique length factorizations if their order class number is 1 or 2.")
    
    non_maximal_HFD_d_values = []
    
    # Combine d values from h=1 and h=2 lists
    d_h1_h2 = d_h1 + d_h2
    
    for d in d_h1_h2:
        # Condition for being a non-maximal order Z[sqrt(-d)]
        if d % 4 == 3:
            # Case d=3
            if d == 3:
                # h(Z[sqrt(-3)])=1, so it's an HFD.
                non_maximal_HFD_d_values.append(d)
            # Case d = 3 (mod 8), d > 3. Condition is h(-d) <= 2.
            # All d in d_h1_h2 satisfy h(-d)<=2, so we just need to check the congruence.
            elif d % 8 == 3:
                non_maximal_HFD_d_values.append(d)
            # Case d = 7 (mod 8) gives no HFDs.
            
    count_non_maximal_HFDs = len(non_maximal_HFD_d_values)
    non_maximal_HFD_d_values.sort()
    
    print(f"\nNumber of such rings: {count_non_maximal_HFDs}")
    print(f"Values of d: {non_maximal_HFD_d_values}")
    
    # Step 3: Total Count
    # The two sets of rings are disjoint, so the total is the sum of the counts.
    total_count = count_maximal_HFDs + count_non_maximal_HFDs

    print("\n--- Final Calculation ---")
    print("The total number of rings is the sum of the counts from the two disjoint sets.")
    print(f"Total = (Integrally Closed HFDs) + (Not Integrally Closed HFDs)")
    print(f"{count_maximal_HFDs} + {count_non_maximal_HFDs} = {total_count}")

    # Return the final numerical answer as requested
    return total_count

if __name__ == '__main__':
    final_answer = solve_unique_length_factorizations()
    # The wrapper format for the final answer
    sys.stdout.write(f"\n<<< {final_answer} >>>\n")