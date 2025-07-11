def solve_iupac_name():
    """
    This function provides the IUPAC name for the compound C7H14 based on the provided 13C NMR data.
    The reasoning is as follows:
    1.  The molecular formula C7H14 has a degree of unsaturation of 1.
    2.  13C NMR signals at 145(s) and 112(t) confirm a >C=CH2 group.
    3.  Having 6 signals for 7 carbons implies two carbons are chemically equivalent.
    4.  The combination of signals at 27(d) and a quartet representing two carbons points to an isopropyl group [-CH(CH3)2].
    5.  The remaining signals at 48(t) and a quartet correspond to an ethyl group [-CH2CH3].
    6.  Assembling these pieces [>C=CH2, ethyl, and isopropyl] yields the structure CH2=C(CH2CH3)(CH(CH3)2).
    7.  The IUPAC name for this structure is 2-ethyl-3-methyl-1-butene.
    """
    
    # The final equation is the IUPAC name, which includes numbers.
    iupac_name = "2-ethyl-3-methyl-1-butene"
    
    print("The IUPAC name of the compound is:")
    print(iupac_name)

solve_iupac_name()
<<<2-ethyl-3-methyl-1-butene>>>