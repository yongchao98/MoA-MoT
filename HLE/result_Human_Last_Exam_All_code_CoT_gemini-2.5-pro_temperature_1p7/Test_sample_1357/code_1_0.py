import math

def main():
    """
    Calculates the number of types of stable reductions for a genus 4 curve
    with a good reduction Jacobian.

    This is equivalent to counting the number of non-isomorphic "stable curves of compact type"
    of genus 4. A stable curve is of "compact type" if its dual graph (components as vertices,
    nodes as edges) is a tree.

    Key conditions for a type:
    1. The sum of the genera of its irreducible components must be 4 (g_1 + ... + g_c = 4).
    2. Any component with genus 0 must connect to at least 3 other components (i.e., its vertex
       in the dual graph must have degree >= 3).

    A "type" is a non-isomorphic tree with vertices weighted by component genera that meets these conditions.
    The script enumerates these types based on the possible partitions of the genus 4.
    """
    total_types = 0
    print("Calculating the number of types of stable reductions for a genus 4 curve with good reduction Jacobian.")
    print("This corresponds to counting non-isomorphic weighted trees representing stable curves of compact type.\n")

    # Cases with no genus 0 components
    # ---------------------------------
    # For these cases, any tree structure is valid. The number of types depends on the genera
    # multiset and the number of non-isomorphic trees for that number of components.

    # Case 1: Partition {4}. One component of genus 4.
    # G = {4}. c=1 component.
    # There is only one tree on 1 vertex.
    count_4 = 1
    total_types += count_4
    print(f"Genera partition {{4}}: {count_4} type (a single g=4 component).")

    # Case 2: Partition {3, 1}. Two components.
    # G = {3, 1}. c=2 components.
    # There is only one tree on 2 vertices. Since genera are distinct (3 and 1), there is only one way to label the tree.
    count_3_1 = 1
    total_types += count_3_1
    print(f"Genera partition {{3, 1}}: {count_3_1} type.")

    # Case 3: Partition {2, 2}. Two components.
    # G = {2, 2}. c=2 components.
    # One tree on 2 vertices. Genera are identical, so there is one unique labeling.
    count_2_2 = 1
    total_types += count_2_2
    print(f"Genera partition {{2, 2}}: {count_2_2} type.")

    # Case 4: Partition {2, 1, 1}. Three components.
    # G = {2, 1, 1}. c=3 components.
    # There is only one tree on 3 vertices (a path). Its vertices are not all symmetric (center vs. leaves).
    # Assigning the g=2 component to the center gives a different type (1-2-1) than assigning it to a leaf (2-1-1).
    count_2_1_1 = 2
    total_types += count_2_1_1
    print(f"Genera partition {{2, 1, 1}}: {count_2_1_1} types.")

    # Case 5: Partition {1, 1, 1, 1}. Four components.
    # G = {1, 1, 1, 1}. c=4 components.
    # Genera are all identical. The number of types is the number of non-isomorphic trees on 4 vertices, which is 2 (a path and a star).
    count_1_1_1_1 = 2
    total_types += count_1_1_1_1
    print(f"Genera partition {{1, 1, 1, 1}}: {count_1_1_1_1} types (path and star graphs).")

    print("\nConsidering cases with genus 0 components (g=0 requires vertex degree >= 3).\n")

    # Case 6: Base partition {2, 1, 1}, plus one g=0 component.
    # G = {2, 1, 1, 0}. c=4 components.
    # We need a tree on 4 vertices where one vertex (for g=0) has degree >= 3.
    # The path graph on 4 vertices has max degree 2, so it's not valid.
    # The star graph has a central vertex of degree 3. g=0 must be at the center. The leaves get genera {2, 1, 1}, which gives one unique labeling.
    count_2_1_1_0 = 1
    total_types += count_2_1_1_0
    print(f"Genera partition {{2, 1, 1, 0}}: {count_2_1_1_0} type.")

    # Case 7: Base partition {1, 1, 1, 1}, plus one g=0 component.
    # G = {1, 1, 1, 1, 0}. c=5 components.
    # We need a tree on 5 vertices where the g=0 vertex has degree >= 3.
    # There are 3 trees on 5 vertices.
    # Tree 1 (path): max degree 2. Invalid.
    # Tree 2 (star, K1,4): has a degree 4 vertex. We must assign g=0 here. Valid. (1 type)
    # Tree 3 ("fork"): has a degree 3 vertex. We must assign g=0 here. Valid. (1 type)
    count_1_1_1_1_0 = 2
    total_types += count_1_1_1_1_0
    print(f"Genera partition {{1, 1, 1, 1, 0}}: {count_1_1_1_1_0} types.")

    # Case 8: Base partition {1, 1, 1, 1}, plus two g=0 components.
    # G = {1, 1, 1, 1, 0, 0}. c=6 components.
    # We need a tree on 6 vertices where two vertices (for g=0) have degree >= 3.
    # Analysis of degree sequences shows this forces the g=1 components to be leaves and the two g=0 components to have degree 3.
    # There is only one such tree structure.
    count_1_1_1_1_0_0 = 1
    total_types += count_1_1_1_1_0_0
    print(f"Genera partition {{1, 1, 1, 1, 0, 0}}: {count_1_1_1_1_0_0} type.")
    
    print("\nOther partitions with more g=0 components are not possible under the degree constraint.")

    print("\n--- Final Calculation ---")
    equation_str = (f"{count_4} (from {{4}}) + "
                    f"{count_3_1} (from {{3,1}}) + "
                    f"{count_2_2} (from {{2,2}}) + "
                    f"{count_2_1_1} (from {{2,1,1}}) + "
                    f"{count_1_1_1_1} (from {{1,1,1,1}}) + "
                    f"{count_2_1_1_0} (from {{2,1,1,0}}) + "
                    f"{count_1_1_1_1_0} (from {{1,1,1,1,0}}) + "
                    f"{count_1_1_1_1_0_0} (from {{1,1,1,1,0,0}})")
    print(equation_str)
    
    final_equation_numbers_only = (f"{count_4} + {count_3_1} + {count_2_2} + {count_2_1_1} + {count_1_1_1_1} + "
                                   f"{count_2_1_1_0} + {count_1_1_1_1_0} + {count_1_1_1_1_0_0} = {total_types}")

    print(final_equation_numbers_only)
    print(f"\nTotal number of types is {total_types}.")

if __name__ == "__main__":
    main()
