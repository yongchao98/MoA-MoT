import re

def identify_compound_c():
    """
    This function provides the identity of compound C based on the reaction scheme.
    The reasoning is based on established chemical transformations:
    1. Condensation to form a xanthylium salt (A).
    2. Nucleophilic substitution of a methoxy group by diethylamine (B).
    3. Complete demethylation of all remaining methoxy ethers to phenols (C).
    """
    compound_name = "1-(diethylamino)-3,6,8-trihydroxy-9-(2,4,6-trihydroxyphenyl)xanthylium"
    molecular_formula = "C23H22NO7" # Formula for the cation, charge is +1

    print(f"The final compound C is: {compound_name}")
    print(f"Its molecular formula is: [{molecular_formula}]+")
    print("\nBreaking down the numbers in the molecular formula:")

    # Use regex to find elements and their counts
    # The pattern finds an uppercase letter, optionally followed by a lowercase letter,
    # and then optionally followed by digits.
    elements = re.findall(r'([A-Z][a-z]?)(\d*)', molecular_formula)

    for element, count in elements:
        # If count is an empty string, it means there's 1 atom of that element
        num_atoms = int(count) if count else 1
        print(f"Number of {element} atoms: {num_atoms}")

identify_compound_c()