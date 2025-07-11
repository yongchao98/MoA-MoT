def solve():
    """
    Calculates the largest value k such that for any valid arrangement of k diamonds
    on a 2024x2024 grid, one diamond can be moved to an adjacent cell.
    """
    n = 2024

    # The problem asks for k = m - 1, where m is the minimum number of diamonds
    # in a "stuck" arrangement.
    # A known construction for a minimal stuck arrangement on an n x n grid (for n even)
    # gives a size of m = (n/2)^2 + n - 1.

    # Number of diamonds in the core lattice (A-type cells)
    n_div_2 = n // 2
    core_diamonds = n_div_2 ** 2

    # Number of diamonds added to the boundary to make the arrangement "stuck"
    boundary_diamonds = n

    # The size of the minimal known stuck arrangement, m
    # One corner diamond is double-counted in the boundary groups
    m = core_diamonds + boundary_diamonds - 1

    # The answer k is m - 1
    k = m - 1

    print("The size of the grid is n x n, where n = {}.".format(n))
    print("A minimal 'stuck' arrangement can be constructed with m diamonds.")
    print("The formula for m is (n/2)^2 + n - 1.")
    print("m = ({}/2)^2 + {} - 1".format(n, n))
    print("m = {}^2 + {} - 1".format(n_div_2, n))
    print("m = {} + {} - 1 = {}".format(core_diamonds, n, m))
    print("\nThe problem asks for the largest k such that any arrangement of k diamonds is not stuck.")
    print("This value is k = m - 1.")
    print("k = {} - 1 = {}".format(m, k))

solve()