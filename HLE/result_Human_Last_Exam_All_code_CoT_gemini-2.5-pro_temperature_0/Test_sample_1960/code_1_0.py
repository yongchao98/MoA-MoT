def generate_formula_string(n, atom="bot"):
    """
    Generates the string representation for a formula made of n atoms.
    f(w) = atom * atom * ... * atom (w times).
    Uses '1' for the multiplicative unit when n=0.
    Uses '*' for the tensor product connective.
    """
    if n < 0:
        raise ValueError("Cannot generate formula for a negative number.")
    if n == 0:
        return "1"
    return " * ".join([atom] * n)

def solve_and_print_encoding():
    """
    Defines a sample equipartitioning problem, calculates the parameters,
    and prints the linear logic encoding.
    """
    # --- Problem Definition ---
    # W: The set of numbers to partition.
    # m: The number of partitions.
    W = {1, 2, 3, 6}
    m = 2

    # --- Calculate b ---
    # b: The target sum for each partition.
    # The condition sum(W) = m * b must hold.
    total_sum = sum(W)
    if total_sum % m != 0:
        print(f"The sum of W ({total_sum}) is not divisible by m ({m}).")
        print("No solution is possible for the equipartitioning problem.")
        return
    b = total_sum // m

    print(f"Equipartitioning Problem: EP(W={W}, m={m}, b={b})")
    print("-" * 20)
    print("The encoding in linear logic is as follows:")
    print("Sequent: {f(w) | w in W} |- C")
    print("-" * 20)

    # --- Print f(w) for each w in W ---
    print("Resource Formulas f(w):")
    # The prompt asks to output each number in the final equation.
    # We interpret this as printing the formula for each number in W.
    for w in sorted(list(W)):
        formula = generate_formula_string(w)
        print(f"f({w}) = {formula}")

    print("-" * 20)

    # --- Print the goal formula C ---
    print("Goal Formula C:")
    partition_formula = generate_formula_string(b)
    # C is m copies of the partition formula, tensored together.
    c_formulas = [f"({partition_formula})" for _ in range(m)]
    c_string = " * ".join(c_formulas)
    print(f"C = {c_string}")

# Execute the function
solve_and_print_encoding()