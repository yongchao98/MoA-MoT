import argparse

def formula_tensor(base, n, connective=" @ "):
    """Creates a formula by tensoring 'base' with itself 'n' times."""
    if n < 0:
        raise ValueError("Cannot create a formula with a negative number of elements.")
    if n == 0:
        return "1"
    return connective.join([base] * n)

def f_to_string(w, atom_alpha='a'):
    """Generates the string for the formula f(w)."""
    return formula_tensor(atom_alpha, w)

def c_to_string(m, b, atom_alpha='a', atom_beta='b'):
    """Generates the string for the formula C(m, b)."""
    if m <= 0 or b <= 0:
        raise ValueError("m and b must be positive integers.")

    # (alpha^b -> beta)
    alpha_b = formula_tensor(atom_alpha, b)
    machine_part = f"({alpha_b} -> {atom_beta})"

    # (machine)^m
    lhs = formula_tensor(machine_part, m)

    # beta^m
    rhs = formula_tensor(atom_beta, m)

    return f"({lhs}) -> ({rhs})"


def main():
    """
    Main function to demonstrate the encoding.
    Provides a default example and allows user input via command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Generate the linear logic encoding for an Equipartitioning Problem instance.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-W', '--weights',
        nargs='+',
        type=int,
        default=[3, 5, 4, 3, 5],
        help="A list of integers for the set W."
    )
    parser.add_argument(
        '-m', '--partitions',
        type=int,
        default=2,
        help="The number of partitions (m)."
    )
    # The script will calculate b based on W and m.
    
    args = parser.parse_args()
    
    W = args.weights
    m = args.partitions
    
    total_sum = sum(W)
    if total_sum % m != 0:
        print(f"Error: The sum of weights ({total_sum}) is not divisible by the number of partitions ({m}).")
        print("The Equipartitioning Problem is not well-defined for this input.")
        return
        
    b = total_sum // m

    print("Equipartitioning Problem Instance:")
    print(f"  W = {W}")
    print(f"  m = {m}")
    print(f"  b = {b} (since sum(W) = {total_sum})")
    print("-" * 30)

    print("Linear Logic Encoding:")
    print("Let 'a' and 'b' be two distinct atomic formulas.")
    print("Let '@' denote the tensor product (⊗) and '->' denote linear implication (⊸).")
    print("-" * 30)

    print("Formulas f(w) for w in W:")
    for w in W:
        print(f"  f({w}) = {f_to_string(w)}")
    print("")

    print("Formula C(m, b):")
    # We print each number in the final formula representation.
    # The numbers are w, m, and b. w's are in f(w), m and b are in C.
    print(f"  C({m}, {b}) = {c_to_string(m, b)}")
    print("")

    print("The sequent {f(w) | w in W} |- C(m, b) is derivable if and only if EP(W, m, b) is true.")


if __name__ == "__main__":
    main()
