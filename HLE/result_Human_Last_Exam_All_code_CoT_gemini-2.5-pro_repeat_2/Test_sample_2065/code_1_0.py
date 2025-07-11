import math

def solve_circulons():
    """
    Calculates the number of circulon defect types for a SO(3) gauge theory
    in d+1 dimensions, for d=1 to 6.
    """
    print("This script calculates the number of 'circulons' (circle-like defects) for a SO(3) gauge theory.")
    print("The classification is given by the order of the group π₁(SO(3)) × π_(d-2)(SO(3)) × π_(d-1)(SO(3)).\n")

    # The order of the homotopy groups π_k(SO(3)).
    # Using math.inf for the infinite group Z.
    # We define the order of π_k for k < 0 to be 1 (trivial group).
    pi_orders = {
        0: 1,
        1: 2,        # π₁(SO(3)) = Z₂
        2: 1,        # π₂(SO(3)) = 0
        3: math.inf, # π₃(SO(3)) = Z
        4: 2,        # π₄(SO(3)) = Z₂
        5: 2         # π₅(SO(3)) = Z₂
    }

    def get_pi_order(k):
        """Returns the order of the k-th homotopy group of SO(3)."""
        if k < 0:
            return 1
        return pi_orders.get(k, "Unknown")

    # Loop through dimensions d=1 to d=6
    for d in range(1, 7):
        if d == 1:
            # Special case: A circle cannot be embedded in R^1.
            # There are no non-trivial defects of this shape.
            # The only configuration is the trivial one.
            print("For d=1: A circle defect cannot exist in 1D space. There is only 1 type (the trivial vacuum).")
            continue

        # Indices for the homotopy groups
        k1 = 1
        k2 = d - 2
        k3 = d - 1

        order1 = get_pi_order(k1)
        order2 = get_pi_order(k2)
        order3 = get_pi_order(k3)

        # Format numbers for printing, replacing inf with '∞'
        f_order1 = '∞' if order1 == math.inf else order1
        f_order2 = '∞' if order2 == math.inf else order2
        f_order3 = '∞' if order3 == math.inf else order3

        # Calculate the total number of circulons
        if any(o == math.inf for o in [order1, order2, order3]):
            result = "infinite"
        else:
            result = order1 * order2 * order3

        # Print the detailed equation and result
        print(f"For d={d}, the number of circulons is |π₁| * |π_{k2}| * |π_{k3}| = {f_order1} * {f_order2} * {f_order3} = {result}.")

solve_circulons()