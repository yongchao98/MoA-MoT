def print_hamiltonicity_threshold_equation():
    """
    This function explains and prints the d-threshold for Hamiltonicity
    for graphs H_n U G(n, p) with minimum degree d >= n/2 - eta.
    """
    print("The d-threshold for Hamiltonicity for a graph H_n U G(n, p)")
    print("with minimum degree delta(H_n) >= n/2 - eta is given by the relation:")
    print()

    # The equation is p = C * eta / n, which we can write as p = C * eta^a * n^b
    equation_str = "p = C * (eta / n)"
    eta_exponent = 1
    n_exponent = 1  # Exponent in the denominator

    print(f"  {equation_str}")
    print()
    print("where C is a sufficiently large positive constant.")
    print("The probability 'p' is proportional to 'eta' and inversely proportional to 'n'.")
    print("\nTo express this with exponents, the equation is p = C * eta^1 * n^(-1).")
    print("The numbers in this final equation are the exponents:")
    print(f"Exponent of eta: {eta_exponent}")
    print(f"Exponent of n: {-n_exponent}")


if __name__ == "__main__":
    print_hamiltonicity_threshold_equation()