import math

def main():
    """
    Calculates the number of types of stable reductions of genus 4 curves
    with good Jacobian reduction.

    The method is a combinatorial enumeration based on the properties of
    stable curves of compact type. The total number of types is the sum of
    the types possible for each number of components (v) in the reduction.
    """
    
    print("This script calculates the number of types of stable reductions for a genus 4 curve with good Jacobian reduction.")
    print("The calculation proceeds by analyzing the possible number of components (v) in the stable reduction.\n")

    # List to store the number of types for each case
    num_types_per_case = []

    # --- Case v=1: 1 component ---
    # Equation: g + d = 4. g is geometric genus, d is number of self-nodes.
    # The pairs (g,d) are (4,0), (3,1), (2,2), (1,3), (0,4).
    # All are stable.
    v1_types = 5
    num_types_per_case.append(v1_types)
    print(f"Number of types with v=1 component: {v1_types}")

    # --- Case v=2: 2 components ---
    # Partition of 4 into 2 parts >= 1: {3,1} and {2,2}.
    # Partition {3,1}:
    #   p=3 types: 4 ((3,0),(2,1),(1,2),(0,3))
    #   p=1 types: 2 ((1,0),(0,1))
    v2_part1_types = 4 * 2
    # Partition {2,2}:
    #   p=2 types: 3 ((2,0),(1,1),(0,2))
    #   Multiset of size 2 from 3 types: (3+2-1 choose 2) = 6
    v2_part2_types = 6
    v2_types = v2_part1_types + v2_part2_types
    num_types_per_case.append(v2_types)
    print(f"Number of types with v=2 components: {v2_types}")

    # --- Case v=3: 3 components ---
    # Partition of 4 into 3 parts >= 1: {2,1,1}.
    # Arrangement 1: (p=1)--(p=2)--(p=1) on a path graph.
    #   Center (p=2): 3 types. Leaves (p=1): multiset of 2 from 2 types -> 3.
    v3_arr1_types = 3 * 3
    # Arrangement 2: (p=2)--(p=1)--(p=1) on a path graph.
    #   Leaf (p=2): 3 types. Center (p=1): 2 types. End leaf (p=1): 2 types.
    v3_arr2_types = 3 * 2 * 2
    v3_types = v3_arr1_types + v3_arr2_types
    num_types_per_case.append(v3_types)
    print(f"Number of types with v=3 components: {v3_types}")
    
    # --- Case v=4: 4 components ---
    # Partition {2,1,1,0}: Realizable on a star graph K_1,3. p=0 is center.
    #   Leaves are {p=2, p=1, p=1}. p=2 leaf: 3 types. p=1 leaves: multiset of 2 from 2 -> 3.
    v4_part1_types = 3 * 3
    # Partition {1,1,1,1}:
    #   On a path graph P_4: 7 distinct configurations (strings of 4 items with reversal symmetry).
    v4_p4_types = 7
    #   On a star graph K_1,3: Center (p=1): 2 types. Leaves (p=1): multiset of 3 from 2 -> 4.
    v4_k13_types = 2 * 4
    v4_types = v4_part1_types + v4_p4_types + v4_k13_types
    num_types_per_case.append(v4_types)
    print(f"Number of types with v=4 components: {v4_types}")

    # --- Case v=5: 5 components ---
    # Partition {1,1,1,1,0}. p=0 component needs degree > 2.
    #   Tree K_1,4 (star graph): p=0 is center. 4 leaves are p=1. Multiset of 4 from 2 -> 5 types.
    v5_k14_types = 5
    #   Tree with degrees (3,2,1,1,1): p=0 is deg=3 vertex. Other 4 are p=1. 12 distinct assignments.
    v5_t3_types = 12
    v5_types = v5_k14_types + v5_t3_types
    num_types_per_case.append(v5_types)
    print(f"Number of types with v=5 components: {v5_types}")

    # --- Case v=6: 6 components ---
    # Partition {1,1,1,1,0,0}. Two p=0 components, need two vertices with degree > 2.
    # Realizable on a tree with degrees (3,3,1,1,1,1).
    # 4 leaves (p=1) attached to the 2 central vertices (p=0).
    # Counting configurations of the 4 leaves gives 6 types.
    v6_types = 6
    num_types_per_case.append(v6_types)
    print(f"Number of types with v=6 components: {v6_types}")
    
    # --- Case v>=7: 0 components ---
    # Partitions of 4 into 7 or more parts must contain at least three 0s.
    # This would require a tree on v>=7 vertices with at least 3 nodes of degree > 2.
    # Such trees do not exist, so there are no types for v>=7.
    v7_types = 0
    num_types_per_case.append(v7_types)
    print(f"Number of types with v>=7 components: {v7_types}")
    
    # --- Final Calculation ---
    total_types = sum(num_types_per_case)
    
    # Filter out the zero case for the final equation string
    equation_parts = [str(n) for n in num_types_per_case if n > 0]
    equation_str = " + ".join(equation_parts)
    
    print("\n-------------------------------------------")
    print("The total number of types is the sum of the types from each case:")
    print(f"{equation_str} = {total_types}")
    print("-------------------------------------------")

if __name__ == "__main__":
    main()
<<<87>>>