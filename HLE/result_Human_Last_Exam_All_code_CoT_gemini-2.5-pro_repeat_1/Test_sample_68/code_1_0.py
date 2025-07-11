import itertools

def solve_knot_coloring():
    """
    This function finds the smallest number of colors needed for a non-trivial
    coloring of the figure-eight knot and demonstrates a valid coloring.

    The figure-eight knot has 4 arcs and 4 crossings. We can label the arcs
    x1, x2, x3, and x4. The coloring rules at each crossing can be represented
    by a system of linear congruences. For the standard figure-eight knot diagram,
    these equations are:
    1. 2*x1 - x2 - x4 = 0 (mod n)
    2. 2*x3 - x1 - x4 = 0 (mod n)
    3. 2*x4 - x2 - x3 = 0 (mod n)
    4. 2*x2 - x1 - x3 = 0 (mod n)

    A knot is n-colorable if its determinant divides n. The determinant of the
    figure-eight knot is 5. Therefore, the smallest n > 1 for a non-trivial
    coloring is 5.

    This script will search for a valid 5-coloring to prove this.
    """
    n = 5
    num_arcs = 4
    colors = range(n)

    print(f"Searching for a non-trivial {n}-coloring of the figure-eight knot...")

    # Generate all possible color assignments for the 4 arcs
    for coloring in itertools.product(colors, repeat=num_arcs):
        x1, x2, x3, x4 = coloring

        # A trivial coloring uses only one color. We are looking for a non-trivial one.
        if len(set(coloring)) > 1:
            # Check if the coloring rules are satisfied at all 4 crossings
            eq1 = (2 * x1 - x2 - x4) % n == 0
            eq2 = (2 * x3 - x1 - x4) % n == 0
            eq3 = (2 * x4 - x2 - x3) % n == 0
            eq4 = (2 * x2 - x1 - x3) % n == 0

            if eq1 and eq2 and eq3 and eq4:
                print("\nFound a valid non-trivial 5-coloring!")
                print(f"The colors for the four arcs are: (x1, x2, x3, x4) = {coloring}")

                print("\nVerifying the coloring equations (mod 5):")

                # Crossing 1: Over-arc x1, Under-arcs x2, x4
                print(f"Crossing 1: 2 * x1 = x2 + x4  =>  2 * {x1} = {x2} + {x4}  =>  {2*x1} = {x2+x4}  =>  {(2*x1)%n} = {(x2+x4)%n}")

                # Crossing 2: Over-arc x3, Under-arcs x1, x4
                print(f"Crossing 2: 2 * x3 = x1 + x4  =>  2 * {x3} = {x1} + {x4}  =>  {2*x3} = {x1+x4}  =>  {(2*x3)%n} = {(x1+x4)%n}")

                # Crossing 3: Over-arc x4, Under-arcs x2, x3
                print(f"Crossing 3: 2 * x4 = x2 + x3  =>  2 * {x4} = {x2} + {x3}  =>  {2*x4} = {x2+x3}  =>  {(2*x4)%n} = {(x2+x3)%n}")

                # Crossing 4: Over-arc x2, Under-arcs x1, x3
                print(f"Crossing 4: 2 * x2 = x1 + x3  =>  2 * {x2} = {x1} + {x3}  =>  {2*x2} = {x1+x3}  =>  {(2*x2)%n} = {(x1+x3)%n}")

                print("\nThe smallest algebraic structure that allows this coloring is the set of colors {0, 1, 2, 3, 4}.")
                print("The number of elements in this structure is 5.")
                return # Stop after finding the first valid coloring

    print("No non-trivial coloring found for the given n.")


if __name__ == '__main__':
    solve_knot_coloring()

<<<5>>>