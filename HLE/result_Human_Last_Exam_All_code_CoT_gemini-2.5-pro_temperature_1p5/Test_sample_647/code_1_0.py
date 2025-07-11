import math

def solve():
    """
    Solves the committee election problem by calculating the satisfaction values
    for W_1 (core) and W_2 (EJR) and then finding their ratio.
    """

    # --- Step 1: Analyze W_1 (Core) and calculate s(N, W_1) ---
    print("### Step 1: Finding s(N, W_1)")
    print("W_1 is a committee of size 20 in the core that minimizes satisfaction for group N = {1..8}.")
    print("\nFirst, we analyze the 'core' property.")
    print("A committee W is in the core if no group of voters can unilaterally deviate to a better committee of the same size, chosen from their commonly approved candidates.")
    print("In this problem, the largest set of commonly approved candidates belongs to group N, with |intersection(A(i) for i in N)| = 8.")
    print("Since the committee size k=20 is larger than 8, no group of voters can form a new committee of size 20.")
    print("Therefore, the core consists of all possible committees of size 20.")
    
    print("\nNext, we find the committee that minimizes s(N, W).")
    print("The satisfaction for group N is s(N, W) = sum_{i in N} |A(i) intersect W|.")
    print("This can be expanded: s(N, W) = 8 * |C_N intersect W| + sum_{i in N} |D_i intersect W|, where C_N = {1..8} and D_i are the two other candidates for voter i in N.")
    print("To minimize this, we should prioritize candidates not approved by any member of N.")
    print("These are the 8 candidates in A(9) U A(10). Let's add them to W_1.")
    print("We need 12 more candidates for our committee of size 20.")
    print("A candidate from C_N adds 8 to the total satisfaction, while a candidate from any D_i adds only 1.")
    print("To keep satisfaction low, we select the remaining 12 candidates from the 16 candidates available in D_1 U D_2 U ... U D_8.")
    
    print("\nWith this construction, W_1 contains 0 candidates from C_N.")
    print("The 12 candidates from the D_i sets each add 1 to the total satisfaction.")
    s_N_W1 = 8 * 0 + 12
    print(f"So, s(N, W_1) = 8 * 0 + 12 = {s_N_W1}")
    
    print("\n" + "="*50 + "\n")

    # --- Step 2: Analyze W_2 (EJR) and calculate s(N, W_2) ---
    print("### Step 2: Finding s(N, W_2)")
    print("W_2 is a committee of size 20 that satisfies Extended Justified Representation (EJR) and minimizes satisfaction for group N.")
    
    print("\nFirst, we analyze the EJR property.")
    print("EJR states that for any l-cohesive group of voters (a group with at least l common candidates), at least one member must approve at least l candidates in the winning committee.")
    print("1. The group {9} is 4-cohesive. EJR implies |A(9) intersect W_2| >= 4. Since |A(9)|=4, W_2 must contain all of A(9).")
    print("2. The group {10} is 4-cohesive. Similarly, W_2 must contain all of A(10).")
    print("   This accounts for 8 candidates in W_2.")
    print("3. The group N={1..8} is 8-cohesive. EJR implies there is a voter i in N with |A(i) intersect W_2| >= 8.")

    print("\nNext, we build W_2 to minimize s(N, W) while satisfying these EJR constraints.")
    print("W_2 contains A(9) U A(10). The other 12 candidates must be chosen from C_N U (U D_i).")
    print("Let m = |C_N intersect W_2|. The other 12-m candidates are from the D_i sets.")
    print("The satisfaction is s(N, W_2) = 8 * m + (12-m) = 7 * m + 12.")
    print("To minimize satisfaction, we must find the minimum possible value for m.")
    
    print("\nThe EJR constraint for group N requires: |A(i) intersect W_2| >= 8 for some i.")
    print("|A(i) intersect W_2| = |C_N intersect W_2| + |D_i intersect W_2| = m + |D_i intersect W_2|.")
    print("Since |D_i|=2, the maximum value for |D_i intersect W_2| is 2.")
    print("So, the EJR condition becomes m + 2 >= 8, which means m must be at least 6.")
    
    print("\nThe minimum possible value for m is 6. This can be achieved by picking 6 candidates for W_2 from C_N, and the remaining 6 from the D_i sets (e.g., D_1, D_2, D_3).")
    m = 6
    s_N_W2 = 7 * m + 12
    print(f"Using the minimal m=6, we get the lowest satisfaction for an EJR committee:")
    print(f"s(N, W_2) = 7 * 6 + 12 = {s_N_W2}")

    print("\n" + "="*50 + "\n")

    # --- Step 3: Calculate the Final Ratio ---
    print("### Step 3: Calculating the Final Ratio")
    ratio = s_N_W1 / s_N_W2
    gcd = math.gcd(s_N_W1, s_N_W2)
    numerator = s_N_W1 // gcd
    denominator = s_N_W2 // gcd
    
    print(f"We have s(N, W_1) = {s_N_W1}")
    print(f"We have s(N, W_2) = {s_N_W2}")
    print("The required ratio is s(N, W_1) / s(N, W_2).")
    print(f"Final Equation: {s_N_W1} / {s_N_W2} = {numerator}/{denominator}")
    
    # Return the final numeric answer for the <<<>>> format
    return ratio

final_ratio = solve()
print(f"<<<{final_ratio}>>>")