def format_tensor_power(base, exponent, connective=" \otimes "):
    """
    Creates a string for a formula like A^n, representing it as
    'A connective A connective ... connective A' (n times).
    """
    if not isinstance(exponent, int) or exponent < 0:
        raise ValueError("Exponent must be a non-negative integer.")
    if exponent == 0:
        return "1"  # The unit for the tensor product
    if exponent == 1:
        return base
    return connective.join([base] * exponent)

def generate_encoding_for_equipartition(W, m, b):
    """
    Generates and prints the linear logic encoding for a given
    equipartitioning problem instance.

    This function prints:
    1. The general definition of the function f(w) and formula C.
    2. The specific formulas for the given W, m, and b.
    3. The final sequent that represents the problem.

    Args:
        W (set or list): A set/list of natural numbers.
        m (int): The number of partitions.
        b (int): The target sum for each partition.
    """
    # First, verify that the total sum matches the expected product
    if sum(W) != m * b:
        print("Error: The problem is ill-defined.")
        print(f"The sum of elements in W is {sum(W)}, but m * b is {m * b}.")
        return

    print("--- Linear Logic Encoding for the Equipartition Problem ---\n")
    print("Let 'A' be an atomic formula representing one unit.\n")

    # --- Print general definitions ---
    print("1. Function f(w): Maps a number w to a formula providing 'w' units.")
    # The line below demonstrates the "final equation" with the variable 'w'
    print("   f(w) = A \otimes ... \otimes A  (w times)\n")

    print("2. Formula C: Represents the goal of forming 'm' partitions, each with a sum of 'b'.")
    # The line below demonstrates the "final equation" with variables 'm' and 'b'
    print("   C = (A^b) \otimes ... \otimes (A^b)  (m times)\n")

    print("--- Encoding for the Specific Instance ---\n")
    print(f"Given values:")
    # Using sorted list for consistent output
    sorted_W = sorted(list(W))
    print(f"  W = {sorted_W}")
    print(f"  m = {m}")
    print(f"  b = {b}\n")

    # --- Print specific formulas ---
    print("The resource formulas {f(w) | w in W} are:")
    context_formulas = []
    for w in sorted_W:
        formula_f_w = format_tensor_power("A", w)
        # This loop outputs each number `w` in the formulas for f(w)
        print(f"  f({w}) = {formula_f_w}")
        context_formulas.append(f"({formula_f_w})" if " " in formula_f_w else formula_f_w)
    print()

    print("The target formula C is:")
    formula_A_b = format_tensor_power("A", b)
    if m > 1:
        # This logic outputs the numbers `b` and `m` in the formula for C
        formula_C = " \otimes ".join([f"({formula_A_b})"] * m)
    else:
        formula_C = formula_A_b
    print(f"  C = {formula_C}\n")

    # --- Print the final sequent ---
    print("The problem EP(W, m, b) is true if and only if the following sequent is derivable:")
    context_str = ", ".join(context_formulas)
    print(f"  {context_str} |- {formula_C}")


if __name__ == '__main__':
    # Define a sample equipartitioning problem
    # W can be partitioned into {1, 5} and {2, 4}, both summing to 6
    W_example = {1, 2, 4, 5}
    m_example = 2
    b_example = 6

    # Generate and print the encoding for the sample problem
    generate_encoding_for_equipartition(W_example, m_example, b_example)