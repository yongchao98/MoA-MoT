def calculate_Dg_sequence():
    """
    Calculates the first 4 terms of the sequence D_g.

    Let g be an integer. Let A_g be the fine moduli stack over the complex
    numbers of principally polarised abelian varieties (PPAVs) of dimension g.
    Let U_g -> A_g be the universal PPAV. D_g is the smallest degree of a
    finite etale cover of A_g such that the polarisation is represented by a
    symmetric line bundle on the base-change of U_g to the cover.

    The value of D_g is the number of even theta characteristics, which is given
    by the formula D_g = 2^(g-1) * (2^g + 1).
    """
    print("The first 4 terms of the sequence D_g are calculated below.")
    print("The formula is: D_g = 2^(g-1) * (2^g + 1)")
    print("-" * 30)

    results = []
    for g in range(1, 5):
        term1 = 2**(g - 1)
        term2 = 2**g + 1
        dg = term1 * term2
        results.append(dg)
        print(f"For g = {g}:")
        print(f"D_{g} = 2^({g}-1) * (2^{g} + 1) = {term1} * {term2} = {dg}")

    print("-" * 30)
    print("The sequence is:", ", ".join(map(str, results)))

if __name__ == "__main__":
    calculate_Dg_sequence()
