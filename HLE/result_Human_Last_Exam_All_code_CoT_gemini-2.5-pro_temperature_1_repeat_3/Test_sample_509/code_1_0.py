def check_section_conditions_for_closed_surface(genus, k):
    """
    Checks the conditions for the existence of a section and a homotopy section
    for the fibration pi_{k, k+1} over a closed, orientable surface M of a given genus.
    
    This illustrates the richer conditions that apply to closed manifolds, as opposed to
    the open manifolds mentioned in the problem, where a section always exists.

    Args:
        genus (int): The genus of the surface M. Must be non-negative.
        k (int): The number of points in the base configuration space. Must be positive.
    """
    if not isinstance(genus, int) or genus < 0:
        print("Error: Genus must be a non-negative integer.")
        return
    if not isinstance(k, int) or k <= 0:
        print("Error: k must be a positive integer.")
        return
        
    print(f"Analyzing fibration for a closed surface of genus g = {genus} and base space conf_{k}(M):")

    # Condition for section existence (Fadell-Neuwirth)
    # A section exists if and only if the Euler characteristic chi(M) is 0.
    # For a closed, orientable surface, chi(M) = 2 - 2g.
    euler_char = 2 - 2 * genus
    print(f"\n--- Section Existence ---")
    print(f"Condition: Euler characteristic chi(M) = 0.")
    print(f"Calculation: chi(M) = 2 - 2*g = 2 - 2*{genus} = {euler_char}")
    if euler_char == 0:
        print(f"Result: Condition met. A section exists.")
        section_exists = True
    else:
        print(f"Result: Condition NOT met. A section does not exist.")
        section_exists = False

    # Condition for homotopy section existence (Bödigheimer-Cohen-Taylor for g>=1, other results for g=0)
    # A homotopy section exists if and only if k >= 2g - 1.
    print(f"\n--- Homotopy Section Existence ---")
    print(f"Condition: k >= 2g - 1")
    rhs = 2 * genus - 1
    print(f"Calculation: Is {k} >= 2*{genus} - 1?  (i.e., Is {k} >= {rhs}?)")
    homotopy_section_exists = (k >= rhs)
    
    if homotopy_section_exists:
        print(f"Result: Condition met. A homotopy section exists.")
    else:
        print(f"Result: Condition NOT met. A homotopy section does not exist.")


# Example 1: A genus 2 surface, which is not an interior of a bounded manifold.
# We check the fibration pi_{2,3}, so k=2.
print("--- Example 1: Genus 2 surface, k=2 ---")
check_section_conditions_for_closed_surface(genus=2, k=2)

# Example 2: The torus (genus 1), for which chi=0.
# We check the fibration pi_{1,2}, so k=1.
print("\n" + "="*40)
print("--- Example 2: The torus (genus 1), k=1 ---")
check_section_conditions_for_closed_surface(genus=1, k=1)

# Example 3: The sphere (genus 0), for which chi != 0.
# We check the fibration pi_{1,2}, so k=1.
print("\n" + "="*40)
print("--- Example 3: The sphere (genus 0), k=1 ---")
check_section_conditions_for_closed_surface(genus=0, k=1)
