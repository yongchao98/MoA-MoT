def generate_and_print_sequent(W, m, b):
    """
    Generates and prints the linear logic sequent that encodes the
    equipartitioning problem EP(W, m, b).

    In the output, '@' represents the tensor connective (⊗) and
    '-o' represents linear implication (⊸).
    """

    # Check if the necessary condition holds
    if sum(W) != m * b:
        print(f"Error: The sum of elements in W ({sum(W)}) is not equal to m * b ({m * b}).")
        print("An equipartition is impossible.")
        return

    def tensor_power(base_formula, n):
        """Constructs a formula by tensoring a base formula n times."""
        if n < 0:
            raise ValueError("Power must be non-negative.")
        if n == 0:
            return "1"
        if n == 1:
            return base_formula
        # Parentheses are important for complex formulas
        # but here we join simple Us, so it's fine.
        return " @ ".join([base_formula] * n)

    # 1. Define the building blocks of the encoding.
    # U is a base formula without any literals, as required.
    U = "(1 -o bot)"

    print("--- Equipartitioning Problem to Linear Logic ---")
    print(f"Problem: Partition W={W} into m={m} subsets each summing to b={b}.")
    print(f"Condition: sum(W) = {sum(W)}, m*b = {m*b}. Condition holds.\n")
    print("--- The Encoding ---\n")
    
    # 2. Define the function f(w) that represents each number.
    print("1. Function f(w): Represents each number w in W as a resource.")
    print(f"   We define a base unit resource U = {U}.")
    print(f"   Then, f(w) is U tensored w times: f(w) = U^@w.")
    
    f_w_formulas = {}
    for w in sorted(list(W)):
        f_w_formulas[w] = tensor_power(U, w)
    
    print("\n   For the given set W, this gives:")
    for w, formula in f_w_formulas.items():
        print(f"   f({w}) = {formula}")

    # 3. Define the goal formula C.
    print("\n2. Goal Formula C: Represents the desired partitioned state.")
    print(f"   C must represent m={m} groups, each with a sum of b={b}.")
    print(f"   Each group is represented by U^@b.")
    u_b_formula = tensor_power(U, b)
    print(f"   U^@{b} = {u_b_formula}")
    print(f"   C is m={m} of these groups tensored together: C = (U^@b)^@m.")
    c_formula = tensor_power(f"({u_b_formula})", m)
    print(f"\n   The final formula for C is:")
    print(f"   C = {c_formula}\n")

    # 4. State the final sequent.
    print("--- The Final Sequent ---\n")
    print("The sequent {f(w) | w in W} |- C is provable iff the equipartition exists.\n")
    
    lhs_str = ", ".join([f"f({w})" for w in sorted(list(W))])
    print(f"{{{lhs_str}}}")
    print("   |-")
    print(f"{c_formula}")


# Example from the thinking process: W={1, 2, 3, 4}, m=2, b=5.
# This problem has a solution: {1, 4} and {2, 3}.
example_W = {1, 2, 3, 4}
example_m = 2
example_b = 5

generate_and_print_sequent(example_W, example_m, example_b)
