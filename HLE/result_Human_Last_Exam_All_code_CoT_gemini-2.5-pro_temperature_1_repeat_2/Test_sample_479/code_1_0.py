import math

def solve_genus_problem():
    """
    This script explains the solution to find the maximal genus of a surface
    with non-vanishing mean curvature.
    """
    print("Problem: Find the maximal genus of a compact surface in R^3 with non-vanishing mean curvature H.")
    print("-" * 30)

    # Step 1: Analyze the condition H != 0
    print("Step 1: The condition H != 0 on a compact, connected surface implies H is either always positive or always negative.")
    print("We can assume, without loss of generality, that H > 0 everywhere.")
    print("-" * 30)

    # Step 2: Check low genus cases
    print("Step 2: Check if surfaces of low genus can satisfy H > 0.")
    
    # Genus 0
    genus_0 = 0
    is_possible_0 = True
    print(f"Case: Genus g = {genus_0}")
    print(f"A sphere is a surface of genus {genus_0}. It has constant positive mean curvature, so it satisfies the condition.")
    print(f"Result: Genus {genus_0} is possible.")
    print()

    # Genus 1
    genus_1 = 1
    is_possible_1 = True
    print(f"Case: Genus g = {genus_1}")
    print(f"A torus is a surface of genus {genus_1}. For a torus of revolution with major radius R and minor radius r,")
    print("the mean curvature can be made positive everywhere if the condition R > 2*r is met.")
    print(f"Result: Genus {genus_1} is possible.")
    print("-" * 30)

    # Step 3: Consider higher genus cases
    print("Step 3: Consider surfaces with genus g >= 2.")
    threshold_genus = 2
    print(f"A known theorem in differential geometry states that any compact, embedded surface in R^3")
    print(f"with genus g >= {threshold_genus} must contain at least one point where the mean curvature H = 0.")
    print("This is because creating handles requires structures like catenoidal necks where H must approach 0.")
    print(f"Result: Genera greater than or equal to {threshold_genus} are NOT possible.")
    print("-" * 30)
    
    # Step 4: Conclusion
    print("Step 4: Conclusion")
    possible_genera = []
    if is_possible_0:
        possible_genera.append(genus_0)
    if is_possible_1:
        possible_genera.append(genus_1)
        
    maximal_genus = max(possible_genera)
    print(f"The possible genera are {possible_genera}.")
    print(f"The maximal possible genus is therefore {maximal_genus}.")

# Run the explanation
solve_genus_problem()