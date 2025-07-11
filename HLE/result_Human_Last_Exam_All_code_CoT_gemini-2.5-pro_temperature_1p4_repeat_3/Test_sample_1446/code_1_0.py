import sys

def calculate_nu_exponent():
    """
    Calculates the critical exponent nu for a phi-4 theory.

    The value depends on the spatial dimension 'd' and the number of
    order parameter components 'n'. This function provides the mean-field result
    for d>=4 and a first-order epsilon expansion for d<4.
    """
    # We select the canonical example of the 3D Ising model universality class.
    # n: number of components of the order parameter
    n = 1
    # d: spatial dimension
    d = 3

    # The upper critical dimension for phi-4 theory is 4.
    d_c = 4

    print(f"Calculating the critical exponent ν for a φ⁴ theory with n={n} and d={d}.")
    print("----------------------------------------------------------------------")

    if d >= d_c:
        nu = 0.5
        print(f"The dimension d={d} is at or above the upper critical dimension (d_c={d_c}).")
        print("Mean-field theory is exact, yielding:")
        print(f"ν = {nu}")
    else:
        # For d < d_c, we use the epsilon expansion.
        epsilon = d_c - d

        print(f"The dimension d={d} is below the upper critical dimension (d_c={d_c}).")
        print("Using the first-order epsilon expansion: ε = d_c - d.")
        print(f"The value of epsilon is: ε = {d_c} - {d} = {epsilon}\n")

        # The first-order formula is: ν ≈ 1/2 + (n+2)/(4*(n+8)) * ε
        c1 = 0.5
        term2_num = n + 2
        term2_den = 4 * (n + 8)
        term2_val = term2_num / term2_den

        nu_approx = c1 + term2_val * epsilon

        print("The first-order approximation formula for ν is:")
        print("ν ≈ 1/2 + ((n + 2) / (4 * (n + 8))) * ε")
        print("\nSubstituting the values:")
        # Show each number in the equation
        print(f"ν ≈ {c1} + (({n} + {2}) / ({4} * ({n} + {8}))) * {epsilon}")
        print(f"ν ≈ {c1} + ({term2_num} / {term2_den}) * {epsilon}")
        print(f"ν ≈ {c1} + {round(term2_val, 4)} * {epsilon}")
        print(f"ν ≈ {c1} + {round(term2_val * epsilon, 4)}")
        print("----------------------------------------------------------------------")
        print(f"The resulting first-order approximation for ν is: {round(nu_approx, 4)}")
        print("\nNote: This is an approximation. Accepted values from higher-order")
        print(f"RG calculations or simulations for the 3D Ising model are ν ≈ 0.630.")

if __name__ == '__main__':
    calculate_nu_exponent()
