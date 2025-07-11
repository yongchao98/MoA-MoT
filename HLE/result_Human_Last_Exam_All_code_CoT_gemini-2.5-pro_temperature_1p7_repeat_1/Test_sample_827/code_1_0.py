def identify_starting_material():
    """
    Analyzes NMR data and reaction pathway to identify the starting material.
    """
    
    # NMR data provided for the final product
    H_NMR = {
        '8.69': 't, J=5.7 Hz, 1H', '8.24': 's, 1H', '8.11': 's, 1H',
        '7.37-7.22': 'm, 5H', '4.73': 'd, J=6.0 Hz, 2H', '1.70': 's, 9H'
    }
    C_NMR = [156.89, 154.96, 152.80, 139.82, 130.16, 128.82, 127.85, 127.35,
             102.23, 59.79, 43.52, 29.25]

    print("Step 1: Analyzing the final product's NMR data.")
    print("-----------------------------------------------")
    print(f"H-NMR signal at δ {list(H_NMR.keys())[3]} ({H_NMR['7.37-7.22']}), δ {list(H_NMR.keys())[4]} ({H_NMR['4.73']}), and δ {list(H_NMR.keys())[0]} ({H_NMR['8.69']}) indicate a benzylamino group (-NH-CH2-Ph).")
    print(f"H-NMR signal at δ {list(H_NMR.keys())[5]} ({H_NMR['1.70']}) confirms a tert-butyl group from tert-butyl hydrazine.")
    print(f"H-NMR signals at δ {list(H_NMR.keys())[1]} ({H_NMR['8.24']}) and δ {list(H_NMR.keys())[2]} ({H_NMR['8.11']}) suggest two uncoupled protons on a heteroaromatic ring.")
    print("The C-NMR data is consistent with a 4,6-disubstituted pyrimidine core with benzylamino and tert-butylhydrazinyl groups.")
    print("Conclusion: The product is 4-(benzylamino)-6-(2-tert-butylhydrazinyl)pyrimidine.")
    print("\n")
    
    print("Step 2: Deducing the starting material (Compound A).")
    print("-----------------------------------------------------")
    print("The reaction is a sequential double nucleophilic aromatic substitution (SNAr) on a heterocyclic ring.")
    print("The two nucleophiles, tert-butyl hydrazine and benzylamine, replace two leaving groups on the starting material.")
    print("Since the product is a 4,6-disubstituted pyrimidine, Compound A must have leaving groups at positions 4 and 6.")
    print("The most common leaving groups for this type of reaction are halogens.")
    print("\n")
    
    print("Final Deduction:")
    print("The starting material (Compound A) is 4,6-dichloropyrimidine.")

if __name__ == "__main__":
    identify_starting_material()
