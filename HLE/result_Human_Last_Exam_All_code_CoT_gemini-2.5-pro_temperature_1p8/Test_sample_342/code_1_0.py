import math

def calculate_Dg_and_explanation(g):
    """
    Calculates the g-th term of the sequence D_g and provides an explanation.

    Args:
        g: The dimension of the principally polarised abelian varieties.

    Returns:
        A tuple containing the value of D_g and a string explaining the calculation.
    """
    # The condition for the obstruction to vanish is g*(g+1)/2 being even,
    # which is equivalent to g being 0 or 3 modulo 4.
    if g % 4 == 0 or g % 4 == 3:
        # If the obstruction vanishes, no cover is needed. D_g = 1.
        return 1, f"For g={g}, g is congruent to {g%4} mod 4, so D_{g} = 1."
    else:
        # If the obstruction is non-trivial, we need a cover.
        # The degree D_g is the number of Lagrangian subspaces in (Z/2Z)^2g.
        product = 1
        factors_str_list = []
        factors_val_list = []

        for i in range(1, g + 1):
            term = 2**i + 1
            product *= term
            factors_str_list.append(f"(2^{i}+1)")
            factors_val_list.append(str(term))

        explanation_str = f"For g={g}, g is congruent to {g%4} mod 4. The value is D_{g} = "
        explanation_str += " * ".join(factors_str_list)
        # Add the intermediate step only if there are multiple factors
        if g > 1:
            explanation_str += f" = {' * '.join(factors_val_list)}"
        explanation_str += f" = {product}"

        return product, explanation_str

def main():
    """
    Main function to calculate and print the first 4 terms of the sequence D_g.
    """
    print("The first 4 terms of the sequence D_g are computed as follows:")
    sequence = []
    for g_val in range(1, 5):
        Dg, desc = calculate_Dg_and_explanation(g_val)
        sequence.append(Dg)
        print(desc)
    print(f"\nThe sequence is: {', '.join(map(str, sequence))}")

if __name__ == "__main__":
    main()