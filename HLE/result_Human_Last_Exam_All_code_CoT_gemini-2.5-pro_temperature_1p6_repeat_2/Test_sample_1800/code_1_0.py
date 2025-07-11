import sys

def calculate_optimal_ratio():
    """
    This function analyzes findings from scientific literature to determine an
    ideal Ni/Ce ratio for catalysis in WGS and WS reactions.

    The catalytic activity of Ni-Ceria is maximized by achieving high dispersion
    of Ni atoms on the ceria support, which creates a large number of active
    Ni-O-Ce interfaces.

    Literature suggests that a nickel loading of around 10 mol% is often optimal.
    - Below this range: There may be an insufficient number of active Ni sites.
    - Above this range: Nickel particles tend to aggregate (sinter) into larger
      clusters, which reduces the active surface area and can block the ceria support.

    This script calculates the atomic Ni/Ce ratio for a representative optimal
    catalyst with 10 mol% Nickel.
    """

    # Define the optimal nickel percentage based on literature consensus.
    ni_mol_percent = 10.0

    # On a basis of 100 total moles of metal atoms (Ni + Ce):
    ni_moles = ni_mol_percent
    ce_moles = 100.0 - ni_mol_percent

    # Calculate the atomic ratio of Ni to Ce
    # This ratio is a key descriptor of the catalyst's composition.
    if ce_moles == 0:
        print("Error: Cerium content cannot be zero.", file=sys.stderr)
        return

    atomic_ratio = ni_moles / ce_moles

    # --- Outputting the result ---
    print("Analysis of the Optimal Ni/Ce Ratio for Ni-Ceria Catalysts")
    print("-" * 60)
    print(f"Based on scientific literature, a common optimal composition is ~{int(ni_mol_percent)} mol% Nickel.")
    print("\nCalculating the corresponding atomic ratio...")

    print("\nStep 1: Define moles on a basis of 100 total metal moles.")
    print(f"  - Moles of Nickel (Ni) = {ni_moles}")
    print(f"  - Moles of Cerium (Ce) = 100 - {ni_moles} = {ce_moles}")

    print("\nStep 2: Calculate the Ni/Ce atomic ratio using the final equation.")
    print(f"  - Final Equation: Ratio = Moles of Ni / Moles of Ce")
    # As requested, printing each number in the final equation during the calculation display
    print(f"  - Calculation: Ratio = {ni_moles} / {ce_moles}")

    print("\n--- Conclusion ---")
    print(f"The calculated ideal Ni/Ce atomic ratio is approximately {atomic_ratio:.2f}.")
    print("This ratio provides a high density of active sites at the Ni-Ce interface while minimizing Ni particle agglomeration, maximizing catalytic performance for reactions like Water Gas Shift.")

if __name__ == '__main__':
    calculate_optimal_ratio()
