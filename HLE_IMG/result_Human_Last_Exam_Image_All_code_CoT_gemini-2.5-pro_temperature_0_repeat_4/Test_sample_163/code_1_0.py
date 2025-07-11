import textwrap

def identify_reaction_products():
    """
    Identifies and describes the two major products (A and B) from the reaction of
    styrene with tert-butyl peroxybenzoate in the presence of an iron catalyst.
    """

    # Information about the reaction and products
    reaction_summary = """
    The reaction of styrene with tert-butyl peroxybenzoate is a radical addition
    across the alkene's double bond. The peroxide cleaves into a benzoyloxy radical
    and a tert-butoxy radical. The addition of these two fragments to styrene can
    occur in two different orders, leading to two major regioisomeric products, A and B.
    """

    product_A = {
        "name": "1-phenyl-2-(tert-butoxy)ethyl benzoate",
        "structure": "Ph-CH(O-C(=O)-Ph)-CH2-O-tBu",
        "smiles": "c1ccccc1C(=O)OC(c2ccccc2)COC(C)(C)C",
        "formation": "Formed by the initial addition of a tert-butoxy radical to styrene, followed by trapping with a benzoyloxy group."
    }

    product_B = {
        "name": "(2-(benzoyloxy)-1-phenylethyl) tert-butyl ether",
        "alt_name": "2-benzoyloxy-1-tert-butoxy-1-phenylethane",
        "structure": "Ph-CH(O-tBu)-CH2-O-C(=O)-Ph",
        "smiles": "c1ccccc1C(OC(C)(C)C)COC(=O)c2ccccc2",
        "formation": "Formed by the initial addition of a benzoyloxy radical to styrene, followed by trapping with a tert-butoxy group."
    }

    # Print the results
    print("Analysis of the reaction products A and B:")
    print(textwrap.fill(reaction_summary, width=80))
    print("\n" + "="*80)

    print("Product A:")
    print(f"  Name: {product_A['name']}")
    print(f"  Structure: {product_A['structure']}")
    print(f"  SMILES: {product_A['smiles']}")
    print(textwrap.fill(f"  Formation: {product_A['formation']}", width=78, initial_indent='  ', subsequent_indent='  '))
    
    print("-" * 80)

    print("Product B:")
    print(f"  Name: {product_B['name']}")
    print(f"  (Alternative Name): {product_B['alt_name']}")
    print(f"  Structure: {product_B['structure']}")
    print(f"  SMILES: {product_B['smiles']}")
    print(textwrap.fill(f"  Formation: {product_B['formation']}", width=78, initial_indent='  ', subsequent_indent='  '))
    
    print("="*80)

if __name__ == '__main__':
    identify_reaction_products()
