import sys
# This script calculates the weight percent (wt%) of Nickel in a Ni-Ceria catalyst
# for a given optimal molar ratio, based on scientific literature.

def calculate_catalyst_composition():
    """
    Calculates and explains the composition of a Ni-Ceria catalyst.

    The ideal Ni/Ce ratio for catalytic activity (WGS, WS) is not a single number
    but a range. Scientific literature suggests an optimal molar fraction of Ni,
    X_Ni = Ni / (Ni + Ce), is often between 0.1 and 0.3. We will use a
    representative value of 0.2 for this calculation.
    """
    # Step 1: Define constants
    # Atomic mass in grams per mole (g/mol)
    atomic_mass_Ni = 58.69
    atomic_mass_Ce = 140.12

    # A representative optimal molar fraction for Ni from literature
    molar_fraction_Ni = 0.2

    # Step 2: Calculate the corresponding molar fraction for Ce
    molar_fraction_Ce = 1 - molar_fraction_Ni

    print("### Catalyst Composition Calculation ###\n")
    print(f"Based on scientific literature, a representative optimal molar fraction of Ni is {molar_fraction_Ni}.")
    print("This means in a sample, for every 1 mole of metal atoms:")
    print(f"- Moles of Nickel (Ni): {molar_fraction_Ni}")
    print(f"- Moles of Cerium (Ce): {molar_fraction_Ce:.1f}\n")

    # Step 3: Calculate the mass of each component based on their molar amounts
    mass_Ni = molar_fraction_Ni * atomic_mass_Ni
    mass_Ce = molar_fraction_Ce * atomic_mass_Ce

    print("--- Calculating Mass from Moles ---")
    print(f"Equation for Ni mass: {molar_fraction_Ni} mol * {atomic_mass_Ni} g/mol = {mass_Ni:.2f} g")
    print(f"Equation for Ce mass: {molar_fraction_Ce:.1f} mol * {atomic_mass_Ce} g/mol = {mass_Ce:.2f} g\n")

    # Step 4: Calculate the total mass and the weight percent (wt%) of Ni
    total_mass = mass_Ni + mass_Ce
    wt_percent_Ni = (mass_Ni / total_mass) * 100

    print("--- Calculating Final Weight Percent (wt%) of Ni ---")
    print(f"Equation for total mass: {mass_Ni:.2f} g (Ni) + {mass_Ce:.2f} g (Ce) = {total_mass:.2f} g")
    print(f"Equation for Ni wt%: ({mass_Ni:.2f} g / {total_mass:.2f} g) * 100 = {wt_percent_Ni:.2f} wt%\n")

    print("### Conclusion ###")
    print(f"A catalyst with a Ni/(Ni+Ce) molar ratio of {molar_fraction_Ni} corresponds to approximately {wt_percent_Ni:.2f} wt% Ni.")
    print("The ideal ratio is a range, but this value is a well-regarded benchmark for achieving high catalytic performance due to strong metal-support interaction.")

if __name__ == '__main__':
    calculate_catalyst_composition()
    # For the final answer, we refer to the optimal molar fraction range of Ni / (Ni + Ce).
    # The value is a range because of varying experimental conditions and synthesis methods.
    sys.stdout.write("\n<<<0.1-0.3>>>\n")