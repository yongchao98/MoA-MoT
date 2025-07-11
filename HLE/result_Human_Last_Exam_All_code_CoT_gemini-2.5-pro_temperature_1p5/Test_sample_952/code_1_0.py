def solve():
    """
    Calculates the largest value k such that for every valid arrangement of k diamonds
    on a 2024x2024 grid, we can move one diamond to an adjacent cell.

    This value k is m_min - 1, where m_min is the minimum size of a "frozen" arrangement.
    A frozen arrangement is one where no single diamond can be moved to an adjacent empty cell.

    For an n x n grid, where n is a multiple of 4 (like 2024), the minimum size of a
    frozen arrangement is known to be (n/2)^2 + 1.
    """
    n = 2024

    # The minimum size of a frozen arrangement for n x n grid where n is a multiple of 4
    # is m_min = (n/2)^2 + 1.
    side_of_subgrid = n // 2
    m_min = side_of_subgrid**2 + 1
    
    # The largest k for which any arrangement is guaranteed to be movable is m_min - 1.
    k = m_min - 1
    
    print(f"The size of the grid is {n}x{n}.")
    print(f"The minimum number of diamonds in a 'frozen' arrangement is m_min = ({n}/2)^2 + 1 = {side_of_subgrid}^2 + 1 = {m_min}.")
    print(f"The largest number of diamonds 'k' that guarantees at least one can be moved is k = m_min - 1.")
    print(f"So, k = {m_min} - 1 = {k}.")
    print(f"The final calculated value is {k}.")

solve()