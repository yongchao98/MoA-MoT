def get_compound_name():
    """
    Determines the IUPAC name of the compound based on spectral data analysis.
    The analysis concludes the structure is 2-phenylpropan-1-amine.
    """
    
    # The structure was determined by analyzing MS, IR, 1H NMR, 13C NMR, DEPT, and HSQC data.
    # The key pieces of evidence are:
    # 1. MS: Molecular ion M+ at m/z 135 and a base peak at m/z 30, indicating a C9H13N formula with a -CH2-NH2 group.
    # 2. 1H & 13C NMR: Reveal a monosubstituted phenyl ring and a propyl side chain with CH, CH2, and CH3 groups.
    # 3. Overall analysis points to 2-phenylpropan-1-amine.
    
    name = "2-phenylpropan-1-amine"
    
    print(f"The IUPAC name of the compound is: {name}")

get_compound_name()