import math

def solve_representation_percentage():
    """
    Calculates the percentage of irreducible representations among all
    indecomposable representations for u_q(sl_2) where q is a primitive
    third root of unity.
    """

    # Let l be the order of the root of unity.
    # q is a primitive third root of unity, so l = 3.
    l = 3

    # Step 1: Count the number of irreducible modules.
    # For a primitive l-th root of unity, there are 'l' non-isomorphic
    # simple (irreducible) modules. Let's call this N_irr.
    num_irreducible = l

    # Step 2: Count the number of reducible indecomposable modules.
    # These fall into two main categories for the small quantum group u_q(sl_2).
    #
    # a) Projective indecomposable modules:
    #    There are (l-1) non-isomorphic projective indecomposable modules
    #    which are the projective covers of the non-projective simple modules.
    num_projective_indecomposable = l - 1
    #
    # b) Periodic modules (sometimes called tadpole modules):
    #    These form l-1 "tubes" in the Auslander-Reiten quiver.
    #    For the standard classification, there are 2 * (l-1) such modules.
    num_periodic_indecomposable = 2 * (l - 1)

    # Step 3: Calculate the total number of indecomposable modules.
    # This is the sum of irreducible and reducible indecomposable modules.
    total_indecomposable_modules = num_irreducible + num_projective_indecomposable + num_periodic_indecomposable

    # Step 4: Calculate the percentage.
    percentage = (num_irreducible / total_indecomposable_modules) * 100

    # Print the detailed breakdown of the calculation.
    print(f"Analysis for u_q(sl_2) with q being a primitive {l}rd root of unity (l=3):")
    print("-" * 60)

    print("Step 1: Counting the irreducible modules (objects in C).")
    print(f"The number of non-isomorphic irreducible modules is l = {num_irreducible}.")

    print("\nStep 2: Counting all indecomposable modules (objects in C).")
    print("The non-isomorphic indecomposable modules are composed of:")
    print(f"  - Irreducible modules: {num_irreducible}")
    print(f"  - Projective indecomposable (reducible) modules: l-1 = {num_projective_indecomposable}")
    print(f"  - Periodic indecomposable (reducible) modules: 2*(l-1) = {num_periodic_indecomposable}")
    
    print("\nThe final equation for the total number of indecomposable modules is:")
    print(f"Total Modules = (Number of Irreducibles) + (Number of Projectives) + (Number of Periodics)")
    print(f"Total Modules = {num_irreducible} + {num_projective_indecomposable} + {num_periodic_indecomposable} = {total_indecomposable_modules}")

    print("\nStep 3: Calculating the final percentage.")
    print("The final equation for the percentage of irreducible modules is:")
    print(f"Percentage = (Number of Irreducible Modules / Total Number of Indecomposable Modules) * 100")
    print(f"Percentage = ({num_irreducible} / {total_indecomposable_modules}) * 100")
    print(f"Result: {percentage:.2f}%")


solve_representation_percentage()