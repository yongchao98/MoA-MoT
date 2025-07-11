def find_ideal_ni_ce_ratio():
    """
    Calculates and explains the ideal Ni/Ce ratio for catalysts based on scientific literature.
    The ideal ratio maximizes the dispersion of Ni active sites on the Ceria support,
    enhancing catalytic activity for Water Gas Shift (WGS) and water splitting reactions.
    """
    # Define the atomic masses (in g/mol).
    atomic_mass_ni = 58.69
    atomic_mass_ce = 140.12

    # Based on literature, an atomic ratio of 1 part Ni to 8 parts Ce is often found
    # to be optimal. This provides a high density of active sites without significant
    # Ni particle agglomeration.
    optimal_ni_atoms = 1
    optimal_ce_atoms = 8

    # Calculate the total mass for this atomic ratio
    total_mass_ni = optimal_ni_atoms * atomic_mass_ni
    total_mass_ce = optimal_ce_atoms * atomic_mass_ce
    total_catalyst_mass = total_mass_ni + total_mass_ce

    # Calculate the weight percentage (wt%) of Nickel in the catalyst
    ni_weight_percentage = (total_mass_ni / total_catalyst_mass) * 100

    # Print the findings
    print("Based on scientific literature, a highly effective Ni/Ce atomic ratio for catalytic performance is:")
    print(f"Ni / Ce = {optimal_ni_atoms} / {optimal_ce_atoms}")
    print(f"\nThis ratio represents a decimal value of {optimal_ni_atoms / optimal_ce_atoms:.3f}.")
    print(f"This atomic ratio corresponds to a Nickel content of approximately {ni_weight_percentage:.2f} wt%.")
    print("\nThis composition is considered ideal as it promotes high dispersion of Ni nanoparticles")
    print("and strong metal-support interaction, which are crucial for maximizing catalytic activity.")

if __name__ == '__main__':
    find_ideal_ni_ce_ratio()