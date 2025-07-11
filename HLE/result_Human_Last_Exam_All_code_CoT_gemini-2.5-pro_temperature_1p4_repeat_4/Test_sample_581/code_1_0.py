import sys

def solve_cap_set_bound():
    """
    This function explains and calculates lower bounds for the size of cap sets
    in dimension 8, and presents the best-known result.
    """

    print("The Cap Set Problem asks for the maximum size of a set of points in an n-dimensional space over the field of 3 elements (F_3^n), such that no three points form a line.")
    print("This size is denoted as r_3(n).\n")

    print("The exact values of r_3(n) are known for small n:")
    known_r3 = {
        1: 2,
        2: 4,
        3: 9,
        4: 20,
        5: 45,
        6: 112
    }
    for n, size in known_r3.items():
        print(f"r_3({n}) = {size}")

    print("\nFor higher dimensions like n=8, the exact value is unknown. We can find lower bounds using constructions.")
    print("A common method is the product construction: r_3(n+m) >= r_3(n) * r_3(m).")
    print("\nLet's apply this to find lower bounds for r_3(8):")

    # Lower bound from r_3(6) and r_3(2)
    n1, m1 = 6, 2
    bound1 = known_r3[n1] * known_r3[m1]
    print(f"Using the partition 8 = {n1} + {m1}, we get r_3(8) >= r_3({n1}) * r_3({m1}).")
    print(f"The resulting equation is: {known_r3[n1]} * {known_r3[m1]} = {bound1}")

    # Lower bound from r_3(5) and r_3(3)
    n2, m2 = 5, 3
    bound2 = known_r3[n2] * known_r3[m2]
    print(f"\nUsing the partition 8 = {n2} + {m2}, we get r_3(8) >= r_3({n2}) * r_3({m2}).")
    print(f"The resulting equation is: {known_r3[n2]} * {known_r3[m2]} = {bound2}")

    # Lower bound from r_3(4) and r_3(4)
    n3, m3 = 4, 4
    bound3 = known_r3[n3] * known_r3[m3]
    print(f"\nUsing the partition 8 = {n3} + {m3}, we get r_3(8) >= r_3({n3}) * r_3({m3}).")
    print(f"The resulting equation is: {known_r3[n3]} * {known_r3[m3]} = {bound3}")

    best_product_bound = max(bound1, bound2, bound3)
    print(f"\nThe best lower bound from these simple product constructions is {best_product_bound}.")

    print("\nHowever, more sophisticated constructions have yielded even better bounds.")
    
    best_known_bound = 496
    print(f"The best known lower bound for the size of a cap set in dimension 8 is {best_known_bound}.")

if __name__ == "__main__":
    solve_cap_set_bound()