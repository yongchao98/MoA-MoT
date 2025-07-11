def calculate_reagent_volume():
    """
    Calculates the precise volume of n-BuLi solution needed for a selective
    lithium-halogen exchange reaction.

    This demonstrates the practical application of using a more precise amount of n-BuLi
    to solve the problem of forming multiple products.
    """

    # --- Chemical & Reaction Parameters ---

    # Molecular weight of 2-bromo-4-chloro-1-iodobenzene (C6H3BrClI)
    MW_SM = 317.35  # g/mol

    # Mass of starting material to be used (example value)
    mass_sm = 5.0  # g

    # Precisely titrated concentration of n-BuLi solution.
    # Accuracy is crucial to prevent side reactions.
    conc_nBuLi = 1.6  # M (mol/L)

    # Desired equivalents of n-BuLi. This is kept very close to 1 to ensure selectivity.
    eq_nBuLi = 1.02

    # --- Stoichiometric Calculation ---

    # 1. Calculate the moles of the starting material
    moles_sm = mass_sm / MW_SM

    # 2. Calculate the target moles of n-BuLi based on the desired equivalents
    moles_nBuLi = moles_sm * eq_nBuLi

    # 3. Calculate the required volume of the n-BuLi solution in milliliters
    # Volume (mL) = (moles / concentration) * 1000
    volume_nBuLi_mL = (moles_nBuLi / conc_nBuLi) * 1000

    # --- Print Results ---

    print("To solve the problem of multiple products, precise reagent addition is key.")
    print("The final calculation is based on the equation:")
    print("Volume (mL) = (Mass_SM / MW_SM) * Equivalents / Concentration * 1000\n")
    print("Plugging in the numbers:")
    # The prompt requires outputting each number in the final equation.
    print(
        f"Volume (mL) = ( {mass_sm:.2f} g / {MW_SM:.2f} g/mol ) * {eq_nBuLi:.2f} / {conc_nBuLi:.2f} M * 1000"
    )
    print(f"Result: {volume_nBuLi_mL:.2f} mL")
    print("\nBy precisely measuring and adding this volume, you supply only enough n-BuLi for the desired reaction,")
    print("thereby preventing the side reaction that creates the second boron-containing product.")

if __name__ == '__main__':
    calculate_reagent_volume()