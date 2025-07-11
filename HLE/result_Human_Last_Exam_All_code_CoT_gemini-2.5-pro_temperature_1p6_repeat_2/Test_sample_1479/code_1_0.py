def describe_fuzzy_mf_structures():
    """
    This function programmatically describes the dimensional structure of
    Type-1, Type-2, and Type-3 fuzzy membership functions to highlight their differences.
    """
    print("--- Analysis of Fuzzy Membership Function (MF) Dimensional Structures ---\n")

    # --- Type-1 MF ---
    t1_domain_variables = 1
    print(f"Type-1 MF: μ(x)")
    print(f"  - Domain Variables: {t1_domain_variables} (primary variable 'x')")
    print(f"  - Description: Maps a single input 'x' to a single membership value.")
    print(f"  - Visualization: 2-Dimensional plot.\n")

    # --- Type-2 MF ---
    t2_domain_variables = 2
    print(f"Type-2 MF: μ(x, u)")
    print(f"  - Domain Variables: {t2_domain_variables} (primary variable 'x', primary membership 'u')")
    print(f"  - Description: Maps an input 'x' to a Type-1 fuzzy set, creating a 3D structure.")
    print(f"  - This is a function with a two-variable domain.\n")

    # --- Type-3 MF ---
    t3_domain_variables = 3
    print(f"Type-3 MF: μ(x, u, v)")
    print(f"  - Domain Variables: {t3_domain_variables} (primary 'x', primary membership 'u', secondary membership 'v')")
    print(f"  - Description: Maps an input 'x' to a Type-2 fuzzy set, creating a 4D structure.")
    print(f"  - This is a function with a three-variable domain.\n")

    # --- Conclusion on the fundamental difference ---
    print("--- Fundamental Difference (Type-2 vs Type-3) ---")
    print(f"A Type-2 MF has a domain of {t2_domain_variables} variables.")
    print(f"A Type-3 MF has a domain of {t3_domain_variables} variables.")
    print("The core structural change is the expansion from a two-variable domain to a three-variable domain.")

describe_fuzzy_mf_structures()