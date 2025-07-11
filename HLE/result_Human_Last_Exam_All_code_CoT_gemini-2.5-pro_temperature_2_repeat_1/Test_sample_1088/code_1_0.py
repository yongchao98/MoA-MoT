def calculate_primitives_6_311G_star_star(element_symbol):
    """
    Calculates and prints the number of primitive Gaussians for an element
    in the 6-311G** basis set. This is valid for first-row atoms.
    """
    # Check if the element is Hydrogen or a "heavy" atom (for this context)
    if element_symbol == 'H':
        # Hydrogen has no core orbitals.
        # Valence 's' primitives from the '311' part:
        valence_s = 3 + 1 + 1
        # Polarization 'p' primitives from the second '*'
        polarization_p = 1
        
        total = valence_s + polarization_p
        
        print(f"For Hydrogen ({element_symbol}):")
        print(f"The total number of primitive Gaussians is {total}.")
        print("Equation:")
        print(f"{valence_s} (valence s) + {polarization_p} (p polarization) = {total}")

    elif element_symbol in ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']:
        # Core 's' primitives from the '6' part:
        core_s = 6
        # Valence primitives from the '311' part, for both s and p orbitals:
        valence_s = 3 + 1 + 1
        valence_p = 3 + 1 + 1
        # Polarization 'd' primitives from the first '*':
        polarization_d = 1

        total = core_s + valence_s + valence_p + polarization_d
        
        print(f"For Carbon ({element_symbol}):")
        print(f"The total number of primitive Gaussians is {total}.")
        print("Equation:")
        print(f"{core_s} (core s) + {valence_s} (valence s) + {valence_p} (valence p) + {polarization_d} (d polarization) = {total}")

    else:
        print(f"The element '{element_symbol}' is not supported by this example.")

# --- Main execution ---
# Calculate and display the results for representative atoms H and C.
print("The 6-311G** basis set composition depends on the atom.")
print("-" * 50)
calculate_primitives_6_311G_star_star('H')
print("-" * 50)
calculate_primitives_6_311G_star_star('C')