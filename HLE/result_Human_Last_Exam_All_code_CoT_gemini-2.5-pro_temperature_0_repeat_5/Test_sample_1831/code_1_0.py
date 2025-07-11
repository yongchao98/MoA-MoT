import math

def main():
    """
    This problem asks for the number of conjugacy classes of endomorphisms on a set of size 4.
    This is equivalent to finding the number of non-isomorphic functional graphs on 4 vertices.

    Our plan is as follows:
    1. A functional graph is a collection of disjoint components. Each component is a
       connected functional graph.
    2. A connected functional graph consists of a cycle with rooted trees attached to its vertices.
    3. We will first find the number of types of connected functional graphs for sizes 1, 2, 3, and 4.
       Let's call this c(n).
    4. Then, we will consider all partitions of the number 4. For each partition, we calculate
       the number of ways to form a graph with components of those sizes.
    5. The sum over all partitions will give the total number of non-isomorphic functional graphs.
    """

    # Step 1: Define t(n), the number of non-isomorphic rooted trees on n vertices.
    # This is a known sequence (OEIS A000081).
    t = {1: 1, 2: 1, 3: 2, 4: 4}
    print("Step 1: Define the number of non-isomorphic rooted trees t(n):")
    print(f"t(1) = {t[1]}, t(2) = {t[2]}, t(3) = {t[3]}, t(4) = {t[4]}\n")

    # Step 2: Calculate c(n), the number of connected functional graphs on n vertices.
    # c(n) is the sum over possible cycle lengths k (from 1 to n) of c(n,k), where c(n,k)
    # is the number of connected graphs on n vertices with a cycle of length k.
    c = {}

    # c(1): Only possible graph is a 1-cycle (a fixed point).
    c[1] = 1
    print("Step 2: Calculate c(n), the number of connected functional graphs on n vertices.")
    print(f"c(1) = {c[1]}")

    # c(2): Cycle of length 1 (a rooted tree t(2)) or cycle of length 2.
    # c(2,1) = t(2) = 1
    # c(2,2) = 1 (a 2-cycle)
    c[2] = t[2] + 1
    print(f"c(2) = (ways for 1-cycle) + (ways for 2-cycle) = {t[2]} + 1 = {c[2]}")

    # c(3): Cycle of length 1 (t(3)), 2, or 3.
    # c(3,1) = t(3) = 2
    # c(3,2): 2-cycle with one vertex attached. This forms a t(2) and a t(1). 1 way.
    # c(3,3): 3-cycle. 1 way.
    c[3] = t[3] + 1 + 1
    print(f"c(3) = (ways for 1-cycle) + (ways for 2-cycle) + (ways for 3-cycle) = {t[3]} + 1 + 1 = {c[3]}")

    # c(4): Cycle of length 1, 2, 3, or 4.
    # c(4,1) = t(4) = 4
    # c(4,2): 2-cycle with 2 vertices to distribute.
    #   - Both on one root: forms a t(3) and a t(1). t(3) has 2 types -> 2 ways.
    #   - One on each root: forms two t(2)s. t(2) has 1 type -> 1 way.
    #   Total for c(4,2) = 2 + 1 = 3
    c4_k2 = t[3] * t[1] + 1
    # c(4,3): 3-cycle with 1 vertex to distribute. Forms one t(2) and two t(1)s. 1 way.
    c4_k3 = 1
    # c(4,4): 4-cycle. 1 way.
    c4_k4 = 1
    c[4] = t[4] + c4_k2 + c4_k3 + c4_k4
    print(f"c(4) = (ways for 1-cycle) + (ways for 2-cycle) + (ways for 3-cycle) + (ways for 4-cycle) = {t[4]} + {c4_k2} + {c4_k3} + {c4_k4} = {c[4]}\n")

    # Step 3: Sum over partitions of 4.
    # The partitions of 4 are: [4], [3,1], [2,2], [2,1,1], [1,1,1,1].
    print("Step 3: Calculate the number of graphs for each partition of 4.")
    
    # Partition [4]: One component of size 4. Number = c(4)
    g4_p4 = c[4]
    print(f"For partition [4]: {g4_p4} graphs.")

    # Partition [3,1]: One component of size 3, one of size 1. Number = c(3) * c(1)
    g4_p31 = c[3] * c[1]
    print(f"For partition [3,1]: c(3) * c(1) = {c[3]} * {c[1]} = {g4_p31} graphs.")

    # Partition [2,2]: Two components of size 2.
    # We need to choose a multiset of size 2 from the c(2)=2 types of components.
    # Number of ways = (n+k-1) choose k, with n=c(2)=2, k=2. (2+2-1) choose 2 = 3.
    g4_p22 = math.comb(c[2] + 2 - 1, 2)
    print(f"For partition [2,2]: {g4_p22} graphs.")

    # Partition [2,1,1]: One component of size 2, two of size 1.
    # There is only c(1)=1 type of size-1 component. So we just choose one component of size 2.
    # Number = c(2)
    g4_p211 = c[2]
    print(f"For partition [2,1,1]: {g4_p211} graphs.")

    # Partition [1,1,1,1]: Four components of size 1.
    # Only c(1)=1 type. So only 1 way.
    g4_p1111 = 1
    print(f"For partition [1,1,1,1]: {g4_p1111} graph.\n")

    # Step 4: Sum the results.
    total = g4_p4 + g4_p31 + g4_p22 + g4_p211 + g4_p1111
    print("Step 4: The total number of elements is the sum of the counts for each partition.")
    print("The final equation is:")
    print(f"{g4_p4} + {g4_p31} + {g4_p22} + {g4_p211} + {g4_p1111} = {total}")

if __name__ == "__main__":
    main()