def calculate_vogel_bound():
    """
    Calculates an upper bound for the braid index of the three-twist knot
    using Vogel's algorithm.
    """
    
    # Properties of the three-twist knot (5_2)
    knot_name = "Three-twist knot (5_2)"
    # Number of crossings in the standard diagram
    c = 5
    # Genus of the knot
    g = 1

    # According to Vogel's algorithm, an upper bound for the braid index is the
    # number of Seifert circles, 's', from the knot's diagram.
    # The number of Seifert circles 's' can be calculated from the number of
    # crossings 'c' and the genus 'g' of the Seifert surface built from
    # that diagram using the formula: s = c - 2g + 1.
    
    s = c - 2 * g + 1

    # --- Output the results ---
    print("Finding the upper bound for the braid index of the three-twist knot using Vogel's algorithm.")
    print("--------------------------------------------------------------------------")
    print(f"Knot: {knot_name}")
    print(f"Number of crossings in the standard diagram (c): {c}")
    print(f"Genus of the knot (g): {g}")
    print("\nAn upper bound for the braid index is the number of Seifert circles (s).")
    print("The formula is: s = c - 2g + 1")
    print("\nCalculation:")
    print(f"s = {c} - (2 * {g}) + 1")
    s_calc_step = c - (2 * g)
    print(f"s = {s_calc_step} + 1")
    print(f"s = {s}")
    print("--------------------------------------------------------------------------")
    print(f"The calculated upper bound for the braid index is {s}.")

if __name__ == '__main__':
    calculate_vogel_bound()