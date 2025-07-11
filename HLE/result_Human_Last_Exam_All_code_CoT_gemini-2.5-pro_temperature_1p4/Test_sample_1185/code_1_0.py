def count_stable_reductions_genus_2():
    """
    Calculates and explains the number of different types of stable reduction
    for curves of genus 2.

    The classification is based on the combinatorial type of the dual graph of the
    stable curve. A curve is stable if its automorphism group is finite.
    """

    print("The types of stable reduction for a genus 2 curve are classified by the structure of the resulting stable curve.")
    print("We can group them by the number of irreducible components.\n")

    # Number of types with 1 irreducible component
    num_types_1_comp = 3
    print("--- Types with 1 irreducible component ---")
    print("1. A smooth curve of genus 2.")
    print("2. An elliptic curve (genus 1) with one node.")
    print("3. A rational curve (genus 0) with two nodes.")
    print(f"Number of types: {num_types_1_comp}\n")

    # Number of types with 2 irreducible components
    num_types_2_comp = 3
    print("--- Types with 2 irreducible components ---")
    print("4. Two elliptic curves meeting at one point.")
    print("5. An elliptic curve and a nodal rational curve meeting at one point.")
    print("6. Two rational curves meeting at three points.")
    print(f"Number of types: {num_types_2_comp}\n")


    # Number of types with 3 irreducible components
    num_types_3_comp = 2
    print("--- Types with 3 irreducible components ---")
    print("7. A chain of an elliptic curve, a rational curve, and another elliptic curve.")
    print("8. A chain of an elliptic curve, a rational curve, and a nodal rational curve.")
    print(f"Number of types: {num_types_3_comp}\n")

    # Calculate and print the total
    total_types = num_types_1_comp + num_types_2_comp + num_types_3_comp
    
    print("--- Total ---")
    print(f"Total number of types = {num_types_1_comp} + {num_types_2_comp} + {num_types_3_comp} = {total_types}")

if __name__ == '__main__':
    count_stable_reductions_genus_2()