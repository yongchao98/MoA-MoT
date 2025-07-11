import numpy as np

def solve_for_active_dimer(Et, Kd):
    """
    Solves for the concentration of active dimer [D] given a total
    enzyme concentration [Et] (expressed as monomer equivalents) and a
    dissociation constant Kd for the equilibrium 2M <=> D.

    The relationship is given by the quadratic equation:
    4[D]^2 - (4[Et] + Kd)[D] + [Et]^2 = 0
    """
    a = 4
    b = -(4 * Et + Kd)
    c = Et**2

    # Solve the quadratic equation
    roots = np.roots([a, b, c])

    # The physically meaningful root is the one where 0 <= [D] <= Et/2
    physical_root = -1
    for root in roots:
        # Check if the root is real and within physical bounds
        if np.isreal(root) and 0 <= root <= Et / 2:
            physical_root = np.real(root)
            break
            
    print(f"For a total enzyme concentration [Et] of {Et:.1f}:")
    # Using format specifiers to show the equation with signs handled correctly
    print(f"  The equilibrium is described by the quadratic equation: {a:.1f}[D]² + ({b:.1f})[D] + {c:.1f} = 0")
    if physical_root != -1:
        percent_active = (2 * physical_root / Et) * 100
        print(f"  The concentration of active dimer [D] is {physical_root:.3f}")
        print(f"  This means {percent_active:.1f}% of the enzyme is in the active form.")
    else:
        print("  No physically meaningful solution found.")
    print("-" * 50)
    
def troubleshoot_assay():
    """
    Demonstrates why increasing enzyme concentration helps overcome cold-induced
    dissociation for an obligate dimer enzyme.
    """
    # Assume a dissociation constant (Kd) that reflects the chilling on ice.
    # A larger Kd means the equilibrium favors the inactive monomers.
    # Units are arbitrary (e.g., uM).
    kd_value = 10.0

    print("--- Troubleshooting Enzyme Assay ---")
    print(f"Hypothesis: Chilling on ice causes the active dimer to dissociate into inactive monomers.")
    print(f"The equilibrium is: 2 * Monomer <=> 1 * Dimer (Active)")
    print(f"We assume a dissociation constant Kd = {kd_value:.1f} under these conditions.\n")

    # Case 1: Low Enzyme Concentration
    low_enzyme_concentration = 2.0
    solve_for_active_dimer(low_enzyme_concentration, kd_value)

    # Case 2: High Enzyme Concentration (10x higher)
    high_enzyme_concentration = 20.0
    solve_for_active_dimer(high_enzyme_concentration, kd_value)
    
    print("\nConclusion:")
    print("Increasing the total enzyme concentration significantly increases the")
    print("percentage of enzyme in the active dimer form. This will shorten the lag")
    print("phase and help produce a linear initial rate for the assay.")
    print("Therefore, increasing enzyme concentration is the correct troubleshooting step.")


if __name__ == '__main__':
    troubleshoot_assay()
