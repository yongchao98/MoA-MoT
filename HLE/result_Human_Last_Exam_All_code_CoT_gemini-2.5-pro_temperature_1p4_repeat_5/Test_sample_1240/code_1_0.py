import numpy as np

def solve_lattice_problems():
    """
    Solves and explains the three questions about lattice root systems.
    """
    
    # --- Introduction to the method ---
    print("This script explains the solutions to the three lattice theory questions.")
    print("The key idea is that for a d-neighbor N of Z^n, the intersection M = Z^n cap N is a sublattice of Z^n.")
    print("M can be described as M = {x in Z^n | x . w = 0 (mod d)} for a primitive integer vector w (the 'glue vector').")
    print("The root system R2(M) consists of vectors v in M with v.v = 2.")
    print("In Z^n, such vectors are of the form v = e_i +/- e_j, where e_k is the k-th standard basis vector.")
    print("So, a root v = e_i +/- e_j is in R2(M) if and only if v . w = w_i +/- w_j is divisible by d.\n")
    
    # --- Question 1 ---
    print("----------------------------------------------------------------------")
    print("1. Is it true that for a d-neighbor N of Z^12, R2(M) can be of type A_11?")
    print("----------------------------------------------------------------------")
    print("(a) Yes\n")
    print("Explanation:")
    print("We need to find a configuration (n, d, w) that produces an A_11 root system.")
    print("The A_11 root system can be represented by vectors {e_i - e_j} for 1 <= i,j <= 12, in a 12-dimensional space.")
    print("This means we need all vectors of the form e_i - e_j to be in R2(M), and no vectors of the form e_i + e_j.\n")

    print("Let's choose our parameters:")
    n_a = 12
    d_a = 3
    w_a = np.ones(n_a, dtype=int)
    print(f"  n = {n_a}")
    print(f"  d = {d_a}")
    print(f"  w = {tuple(w_a)} (This vector is primitive as gcd is 1)")

    print("\nNow we check the conditions for roots v = e_i +/- e_j:")
    print(" - For v = e_i - e_j (with i != j):")
    print(f"   v . w = w_i - w_j = {w_a[0]} - {w_a[1]} = 0.")
    print(f"   Is 0 divisible by d={d_a}? Yes. So, all e_i - e_j are in R2(M).")
    
    print("\n - For v = e_i + e_j (with i != j):")
    print(f"   v . w = w_i + w_j = {w_a[0]} + {w_a[1]} = 2.")
    print(f"   Is 2 divisible by d={d_a}? No.")
    print("   If we had picked d=2, this would fail. With d=3, it works. So, no e_i + e_j are in R2(M).")

    print("\nConclusion: The resulting root system R2(M) is precisely {e_i - e_j | i,j in {1..12}}, which is the definition of A_11.")

    # --- Question 2 ---
    print("\n--------------------------------------------------------------------------------------------")
    print("2. Can the visible root system R2(M) of a d-neighbor N of Z^15 contain a D_7 component?")
    print("--------------------------------------------------------------------------------------------")
    print("(b) yes\n")
    print("Explanation:")
    print("A D_7 component requires that for a subset of 7 indices, say I = {1, ..., 7}, all vectors {e_i +/- e_j} for i,j in I are in R2(M).")
    print("This means we need w_i + w_j and w_i - w_j to be divisible by d for all i,j in I.\n")

    print("Let's choose our parameters:")
    n_b = 15
    d_b = 2
    w_b = np.zeros(n_b, dtype=int)
    I_b = range(7) # Indices for the D_7 component
    for i in I_b:
      w_b[i] = 1
    
    print(f"  n = {n_b}")
    print(f"  d = {d_b}")
    print(f"  w = {tuple(w_b)} (Primitive, as it contains 1s)")

    print("\nNow we check the conditions for roots v = e_i +/- e_j, where i,j are in {1, ..., 7}:")
    print(" - For v = e_i - e_j:")
    print(f"   v . w = w_i - w_j = {w_b[0]} - {w_b[1]} = 0.")
    print(f"   Is 0 divisible by d={d_b}? Yes.")
    
    print("\n - For v = e_i + e_j:")
    print(f"   v . w = w_i + w_j = {w_b[0]} + {w_b[1]} = 2.")
    print(f"   Is 2 divisible by d={d_b}? Yes (since 2 % 2 == 0).")
    
    print("\nConclusion: All vectors e_i +/- e_j for i,j in {1,...,7} are in R2(M). This forms a D_7 component.")
    print("(Note: R2(M) for this construction also contains a D_8 component on the remaining indices where w_k=0, so R2(M) = D_7 + D_8).")


    # --- Question 3 ---
    print("\n--------------------------------------------------------------------------------------------")
    print("3. For n = 18 and d = 5, is it possible for R2(M) to include more than one D_n component?")
    print("--------------------------------------------------------------------------------------------")
    print("(c) yes\n")
    print("Explanation:")
    print("A D_k component on an index set I requires w_i +/- w_j = 0 (mod d) for all i,j in I.")
    print("This implies 2*w_i = 0 (mod d) and w_i = w_j (mod d).")
    print("For d=5, 2*w_i = 0 (mod 5) implies w_i = 0 (mod 5), since gcd(2,5)=1.")
    print("So, any D_k component must be built on indices 'i' where w_i is a multiple of 5.\n")
    
    print("The root system R2(M) decomposes into orthogonal components based on the values of w_i (mod d).")
    print("For d=5, the structure is: D_{n_0} + A_{n_1+n_4-1} + A_{n_2+n_3-1}, where n_k is the number of indices i with w_i = k (mod 5).")
    print("To get more than one D-type component, we can use the fact that A_3 is isomorphic to D_3.")
    print("We can try to construct a configuration where we have a D_{n_0} component and one of the A-type components becomes a D_3.\n")

    print("Let's engineer the counts n_k to get D_2 and D_3 components:")
    n_k_counts = {'n0': 2, 'n1': 12, 'n2': 2, 'n3': 2, 'n4': 0}
    print(f"  Target counts: n_0=2 (for a D_2), n_2=2 and n_3=2 (to make A_{{2+2-1}} = A_3 ~= D_3). n_1=12 to sum to 18.")
    
    print("\nLet's choose our parameters:")
    n_c = 18
    d_c = 5
    w_c = np.ones(n_c, dtype=int)
    w_c[0:2] = 2  # S_2
    w_c[2:4] = 3  # S_3
    w_c[4:6] = 0  # S_0
    # The rest are 1 (S_1)
    
    print(f"  n = {n_c}")
    print(f"  d = {d_c}")
    print(f"  w = {tuple(w_c)} (Primitive, gcd(0,1,2,3)=1)")

    print("\nAnalyzing the components of the root system R2(M):")
    print(f" - Component on indices where w_i=0 (mod 5) (indices 5,6): This forms a D_{n_0} = D_2 component.")
    print(f" - Component on indices where w_i=2 or 3 (mod 5) (indices 1,2,3,4): This forms an A_{{n_2+n_3-1}} = A_{{2+2-1}} = A_3 component.")
    print(f" - Component on indices where w_i=1 or 4 (mod 5) (indices 7-18): This forms an A_{{n_1+n_4-1}} = A_{{12+0-1}} = A_{11} component.")

    print("\nResulting root system R2(M) = D_2 + A_3 + A_11.")
    print("Since A_3 is isomorphic to D_3, the root system can be written as D_2 + D_3 + A_11.")
    print("Conclusion: This system contains a D_2 and a D_3 component, so it has more than one D_n component.")
    
    # --- Final Answer ---
    print("\n\n<<< (a) Yes; (b) yes; (c) yes >>>")

solve_lattice_problems()