def solve_stable_reduction_genus_2():
    """
    This function determines and explains the number of different types of stable reduction
    for curves of genus 2.
    """
    # The types of stable reduction for a genus 2 curve correspond to the topological
    # types of stable curves of arithmetic genus 2. We can list them based on their
    # components and singularities.

    # 1. Irreducible types (one component)
    irreducible_types = [
        "A smooth irreducible curve of geometric genus 2.",
        "An irreducible curve of geometric genus 1 with one node.",
        "An irreducible curve of geometric genus 0 with two nodes."
    ]
    num_irreducible = len(irreducible_types)

    # 2. Reducible types (more than one component)
    reducible_types = [
        "Two curves of genus 1 intersecting at one point.",
        "Two curves of genus 0 intersecting at three points."
    ]
    num_reducible = len(reducible_types)

    # Calculate the total number of types
    total_types = num_irreducible + num_reducible

    print(f"There are {total_types} different types of stable reduction for curves of genus 2.")
    print("\nThey are classified as follows:")
    
    print("\nIrreducible Types:")
    for i, desc in enumerate(irreducible_types, 1):
        print(f"{i}. {desc}")
        
    print("\nReducible Types:")
    for i, desc in enumerate(reducible_types, 1):
        print(f"{i}. {desc}")

    print("\nThe total number is the sum of the counts for each category.")
    print(f"Final calculation: {num_irreducible} + {num_reducible} = {total_types}")

solve_stable_reduction_genus_2()