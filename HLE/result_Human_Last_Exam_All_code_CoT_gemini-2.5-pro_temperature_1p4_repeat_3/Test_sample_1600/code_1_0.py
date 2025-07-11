def solve_feynman_sequence(target_n):
    """
    Calculates a(n), the number of non-vanishing Feynman diagrams of order 2n,
    using the recurrence a(n) = a(n-1) + (n-1)*a(n-2).
    """
    if not isinstance(target_n, int) or target_n < 0:
        print("Error: Input must be a non-negative integer.")
        return

    # Initialize a list to store the sequence values, starting with a(0) and a(1)
    # a[i] will store the value of a(i)
    a = [0] * (target_n + 1)
    if target_n >= 0:
        a[0] = 1
    if target_n >= 1:
        a[1] = 1

    # Iteratively compute the sequence up to the target number
    for i in range(2, target_n + 1):
        a[i] = a[i-1] + (i-1) * a[i-2]

    # If the target is a(3), display the specific calculation as requested
    if target_n == 3:
        val_a2 = a[2]
        val_a1 = a[1]
        val_a3 = a[3]
        
        print(f"To find a(3), we use the recurrence relation a(n) = a(n-1) + (n-1)*a(n-2).")
        print(f"For n=3, the formula is: a(3) = a(2) + (3-1)*a(3-2) = a(2) + 2*a(1).")
        print(f"We have previously calculated a(1) = {val_a1} and a(2) = {val_a2}.")
        print(f"Substituting these values into the equation for a(3):")
        # Output each number in the final equation
        print(f"{val_a2} + 2 * {val_a1} = {val_a3}")
    else:
        print(f"The value of a({target_n}) is: {a[target_n]}")

# The user asks for a(3)
solve_feynman_sequence(3)