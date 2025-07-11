def solve_cap_set_bound():
    """
    This function calculates a historical lower bound for the cap set problem
    in dimension 8 and then states the current best known bound.
    """
    # The maximum size of a cap set in dimension n is denoted by r_3(n).
    # The exact values for dimensions 3 and 5 are known.
    r3_3 = 9
    r3_5 = 45

    # A lower bound can be found using the product construction:
    # r_3(n1 + n2) >= r_3(n1) * r_3(n2)
    # For n=8, we can use n1=3 and n2=5.
    product_bound = r3_3 * r3_5

    print("A well-known method to find a lower bound for the size of a cap set in dimension 8 is the product construction.")
    print(f"Using the known sizes for dimensions 3 and 5, this method gives the equation:")
    print(f"{r3_3} * {r3_5} = {product_bound}")
    print(f"This provided a lower bound of {product_bound}.")

    print("\nHowever, in 2017, a new construction by Edel and Wassermann established a better lower bound.")

    # The current best known lower bound from the 2017 paper by Edel and Wassermann.
    best_known_bound_8d = 496
    print(f"The best known lower bound for the size of cap sets in dimension 8 is: {best_known_bound_8d}")

solve_cap_set_bound()