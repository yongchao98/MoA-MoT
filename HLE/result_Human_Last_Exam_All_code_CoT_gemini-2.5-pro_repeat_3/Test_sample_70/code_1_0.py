def find_rank_2_flags():
    """
    Identifies African national flags with a linear algebraic rank of 2.

    The rank of a flag's matrix representation is determined by its geometric
    complexity.
    - Rank 1: Simple patterns like solid colors or parallel stripes.
    - Rank 2: A simple pattern combined with a second feature like an emblem,
              canton, or hoist triangle. The Danish flag's cross is rank 2.
    - Rank 3+: Complex patterns with multiple intersecting or non-aligned
               features.

    This function classifies each African flag's structure and filters for
    those with a rank of 2.
    """
    # A dictionary mapping each African country to its flag's structural type.
    african_flags_structure = {
        "Algeria": "vertical_stripes_with_emblem",
        "Angola": "horizontal_stripes_with_emblem",
        "Benin": "vertical_pale_on_horizontal_stripes",
        "Botswana": "horizontal_stripes",
        "Burkina Faso": "horizontal_stripes_with_emblem",
        "Burundi": "saltire_cross",
        "Cabo Verde": "horizontal_stripes_with_emblem",
        "Cameroon": "vertical_stripes_with_emblem",
        "Central African Republic": "stripes_and_pale_and_canton",
        "Chad": "vertical_stripes",
        "Comoros": "horizontal_stripes_with_hoist_triangle",
        "Democratic Republic of the Congo": "complex_diagonal_stripe_with_canton",
        "Republic of the Congo": "complex_diagonal_stripe",
        "Djibouti": "horizontal_stripes_with_hoist_triangle",
        "Egypt": "horizontal_stripes_with_emblem",
        "Equatorial Guinea": "horizontal_stripes_with_hoist_triangle",
        "Eritrea": "complex_hoist_triangle",
        "Eswatini": "horizontal_stripes_with_emblem",
        "Ethiopia": "horizontal_stripes_with_emblem",
        "Gabon": "horizontal_stripes",
        "Gambia": "horizontal_stripes",
        "Ghana": "horizontal_stripes_with_emblem",
        "Guinea": "vertical_stripes",
        "Guinea-Bissau": "vertical_pale_on_horizontal_stripes",
        "Ivory Coast": "vertical_stripes",
        "Kenya": "horizontal_stripes_with_emblem",
        "Lesotho": "horizontal_stripes_with_emblem",
        "Liberia": "horizontal_stripes_with_canton",
        "Libya": "horizontal_stripes_with_emblem",
        "Madagascar": "vertical_pale_on_horizontal_stripes",
        "Malawi": "horizontal_stripes_with_emblem",
        "Mali": "vertical_stripes",
        "Mauritania": "horizontal_stripes_with_emblem",
        "Mauritius": "horizontal_stripes",
        "Morocco": "monocolor_with_emblem",
        "Mozambique": "horizontal_stripes_with_hoist_triangle",
        "Namibia": "complex_diagonal_stripe_with_canton",
        "Niger": "horizontal_stripes_with_emblem",
        "Nigeria": "vertical_stripes",
        "Rwanda": "horizontal_stripes_with_emblem",
        "São Tomé and Príncipe": "horizontal_stripes_with_hoist_triangle",
        "Senegal": "vertical_stripes_with_emblem",
        "Seychelles": "fanned_rays",
        "Sierra Leone": "horizontal_stripes",
        "Somalia": "monocolor_with_emblem",
        "South Africa": "y_shape_pall",
        "South Sudan": "horizontal_stripes_with_hoist_triangle",
        "Sudan": "horizontal_stripes_with_hoist_triangle",
        "Tanzania": "complex_diagonal_stripe",
        "Togo": "horizontal_stripes_with_canton",
        "Tunisia": "monocolor_with_emblem",
        "Uganda": "horizontal_stripes_with_emblem",
        "Zambia": "monocolor_with_fly_canton",
        "Zimbabwe": "horizontal_stripes_with_hoist_triangle",
    }

    # Define which structural categories correspond to a rank of 2.
    # These are structures with two independent geometric features.
    structures_rank_2 = {
        "monocolor_with_emblem",
        "horizontal_stripes_with_emblem",
        "vertical_stripes_with_emblem",
        "horizontal_stripes_with_canton",
        "horizontal_stripes_with_hoist_triangle",
        "vertical_pale_on_horizontal_stripes",
        "monocolor_with_fly_canton",
    }

    # Filter the countries to find those with rank 2 flags.
    rank_2_countries = []
    for country, structure in african_flags_structure.items():
        if structure in structures_rank_2:
            rank_2_countries.append(country)

    # Sort the list for clean presentation and print the result.
    rank_2_countries.sort()
    print("The following African nations have flags with the same linear algebraic rank (Rank 2) as the flag of Denmark:")
    for country in rank_2_countries:
        print(f"- {country}")

find_rank_2_flags()