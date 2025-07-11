import sympy

def identify_manifold():
    """
    This function outlines the steps to identify the 3-manifold
    from the given Heegaard diagram by calculating its fundamental group.
    """
    # Define symbolic generators for the fundamental group
    a_1, a_2, a_3 = sympy.symbols('a_1 a_2 a_3')

    print("Step 1: Identify the Heegaard diagram components.")
    print("The diagram is a genus-3 Heegaard diagram.")
    print("The red curves (α_1, α_2, α_3) define the generators of the fundamental group.")
    print(f"Generators: {a_1}, {a_2}, {a_3}\n")

    print("Step 2: Identify the relators from the blue curves (β_1, β_2, β_3).")
    print("The blue graph is a triangular prism. The β-curves are the boundaries of its three quadrilateral faces.")
    print("We derive the relators by expressing each β-curve as a word in the generators.\n")

    print("Step 3: Derive the relator equations.")
    # Relator 1
    # The numbers in the equation are the indices of the generators.
    print(f"Relator R_1 from the loop around handles 1 and 2 gives the relation: a_2 * a_1**(-1) = 1")
    print(f"  => a_2 = a_1")
    # Relator 2
    print(f"Relator R_2 from the loop around handles 2 and 3 gives the relation: a_3 * a_2**(-1) = 1")
    print(f"  => a_3 = a_2")
    # Relator 3
    print(f"Relator R_3 from the loop around handles 3 and 1 gives the relation: a_1 * a_3**(-1) = 1")
    print(f"  => a_1 = a_3\n")


    print("Step 4: Determine the fundamental group π_1(M).")
    print("The group presentation is <a_1, a_2, a_3 | a_2 = a_1, a_3 = a_2, a_1 = a_3>.")
    print("All three relations are consistent and imply a_1 = a_2 = a_3.")
    print("The group collapses to a single generator, let's call it 'a', with no relations.")
    print("So, π_1(M) = <a>, which is the infinite cyclic group, Z.\n")

    print("Step 5: Identify the 3-manifold.")
    print("The fundamental group π_1(M) = Z is a strong invariant.")
    print("The only prime 3-manifold with this fundamental group is S^2 x S^1 (the product of a 2-sphere and a circle).\n")

    manifold_name = "S^2 x S^1"
    print(f"Conclusion: The Heegaard diagram represents the 3-manifold {manifold_name}.")

if __name__ == "__main__":
    identify_manifold()