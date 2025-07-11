import sys

def find_ideal_ratio():
    """
    This function explains and provides the ideal Ni/Ce ratio for catalysis
    based on findings from scientific literature.
    """
    print("The ideal Ni/Ce ratio for Ni-Ceria catalysts in Water Gas Shift (WGS) and water splitting reactions is determined by a balance between two factors:")
    print("1. Maximizing the number of active sites, which are located at the interface between Ni and Ceria (CeO2).")
    print("2. Preventing the formation of large, bulk Nickel Oxide (NiO) clusters, which are less catalytically active and prone to deactivation (sintering).\n")

    print("Experimental studies have consistently shown that a relatively low concentration of nickel is optimal.")
    print("This ensures that the nickel is highly dispersed as small nanoparticles on the ceria support, maximizing the crucial Ni-CeO2 interface area.\n")

    # Define the optimal range based on literature
    min_ni_atomic_percent = 10
    max_ni_atomic_percent = 20

    # Calculate the corresponding molar ratios
    # Ratio = Ni / Ce = Ni / (100 - Ni)
    min_ratio = min_ni_atomic_percent / (100 - min_ni_atomic_percent)
    max_ratio = max_ni_atomic_percent / (100 - max_ni_atomic_percent)

    print(f"The optimal range for the Nickel (Ni) content is generally found to be between {min_ni_atomic_percent}% and {max_ni_atomic_percent}% of the total metal atoms (Ni + Ce).")
    print(f"This corresponds to a Ni/Ce molar ratio in the range of {min_ratio:.2f} to {max_ratio:.2f}.")
    print("\nA Ni content within this range provides a high density of active sites while maintaining catalyst stability.")

    # The final answer is formatted and written to stdout.
    # The format <<<answer>>> is a special instruction for the system.
    final_answer = f"A Ni atomic percentage of {min_ni_atomic_percent}-{max_ni_atomic_percent}% (a Ni/Ce molar ratio of approximately {min_ratio:.2f} to {max_ratio:.2f})"
    sys.stdout.write(f"\n<<<{final_answer}>>>")

if __name__ == '__main__':
    find_ideal_ratio()