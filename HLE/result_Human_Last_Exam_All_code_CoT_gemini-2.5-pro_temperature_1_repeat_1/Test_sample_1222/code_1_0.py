def check_quiver_taft_condition(n, d):
    """
    Checks the condition on d for a quiver-Taft map σ to potentially be
    non-zero on all arrows.

    The condition is that n-d must be odd.

    Args:
        n (int): The number of vertices, used in the reflection formula.
        d (int): The shift parameter for the reflection.
    """
    # The condition is that n-d is odd.
    value = n - d
    is_odd = (value % 2) != 0

    print(f"For the given parameters n = {n} and d = {d}:")
    print("The condition for σ(a) to be non-zero for all arrows 'a' is that 'n-d' must be odd.")
    print(f"The final equation to check is: ({n} - {d}) % 2 != 0")
    print(f"Calculating the value of n - d: {n} - {d} = {value}")
    print(f"Is {value} odd? {is_odd}")
    if is_odd:
        print("The condition is satisfied.")
    else:
        print("The condition is not satisfied.")

# Example usage:
# Case 1: Condition is satisfied
print("--- Case 1 ---")
check_quiver_taft_condition(n=10, d=3)

print("\n--- Case 2 ---")
# Case 2: Condition is not satisfied
check_quiver_taft_condition(n=10, d=4)