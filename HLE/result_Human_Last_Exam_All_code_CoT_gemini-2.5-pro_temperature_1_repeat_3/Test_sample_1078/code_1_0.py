import math

def main():
    """
    This script identifies and counts the non-isomorphic finite Weyl groups of rank 4.
    """

    # The non-isomorphic finite Weyl groups of rank 4 are determined by the
    # finite-type Dynkin diagrams with 4 vertices.
    # These are A_4, B_4, C_4, D_4, and F_4.
    # The Weyl groups for B_n and C_n are isomorphic, so they count as one type.
    # The other types (A_4, D_4, F_4) are non-isomorphic to each other or to B_4.
    # We can confirm they are non-isomorphic by checking their orders.

    groups = {
        "A_4": {
            "count": 1,
            "order_formula": "5!",
            "order_value": math.factorial(5)
        },
        "B_4 / C_4": {
            "count": 1,
            "order_formula": "2^4 * 4!",
            "order_value": (2**4) * math.factorial(4)
        },
        "D_4": {
            "count": 1,
            "order_formula": "2^(4-1) * 4!",
            "order_value": (2**(4-1)) * math.factorial(4)
        },
        "F_4": {
            "count": 1,
            "order_formula": "1152",
            "order_value": 1152
        }
    }

    print("The non-isomorphic finite Weyl groups of rank 4 and their orders are:")
    for name, data in groups.items():
        print(f"- Type {name}: Order = {data['order_value']}")

    # Create and print the equation for the total count.
    counts = [str(data['count']) for data in groups.values()]
    total_count = sum([data['count'] for data in groups.values()])

    equation = " + ".join(counts) + f" = {total_count}"

    print("\nEach type contributes 1 to the total count of non-isomorphic groups.")
    print("The final calculation is:")
    print(equation)

if __name__ == "__main__":
    main()