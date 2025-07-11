def compare_lipid_surface_area():
    """
    Compares the limiting surface area per molecule for C16-dihydroceramide
    and C16-ceramide based on typical experimental values.

    The area is measured in Angstroms squared (Å²) per molecule.
    """
    # Based on experimental data from pressure-area isotherms, we can assign
    # approximate limiting molecular areas.
    # C16-dihydroceramide (saturated chains) packs tightly.
    # C16-ceramide (one trans-double bond) has a kink and packs less tightly.
    lipid_areas = {
        "C16-dihydroceramide": 40.0,  # Smaller area due to tight packing
        "C16-ceramide": 44.0        # Larger area due to the double bond kink
    }

    # Find the lipid with the minimum area
    lipid_with_lower_area = min(lipid_areas, key=lipid_areas.get)
    lower_area_value = lipid_areas[lipid_with_lower_area]

    lipid_with_higher_area = max(lipid_areas, key=lipid_areas.get)
    higher_area_value = lipid_areas[lipid_with_higher_area]

    print("Comparing molecular surface areas in a compressed monolayer:")
    print(f"- {lipid_with_higher_area} (less ordered) has a surface area of approximately {higher_area_value} Å²/molecule.")
    print(f"- {lipid_with_lower_area} (highly ordered) has a surface area of approximately {lower_area_value} Å²/molecule.")
    print("\nConclusion:")
    print(f"The lipid with the lower surface area is {lipid_with_lower_area}.")

if __name__ == "__main__":
    compare_lipid_surface_area()