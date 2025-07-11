import sys

def calculate_nmr_peaks():
    """
    Calculates the expected number of peaks in the 1H NMR spectrum of
    1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene.

    The calculation is based on molecular symmetry and the number of chemically non-equivalent protons.
    """
    # Step 1: Analyze the molecular symmetry.
    # The molecule has C3 rotational symmetry around the center of the benzene ring.
    # This makes the three -CH3 groups on the ring equivalent, and the three large substituent arms equivalent.

    # Step 2: Count peaks from the central core.
    # The three methyl groups on the central benzene ring are equivalent.
    peaks_central_methyls = 1
    
    # Step 3: Count peaks from one of the three identical substituent arms.
    # Each arm is -CH2-R where R is the chiral indazole moiety.
    
    # a) The -CH2- linker group.
    # Since the attached indazole group (R) is chiral, the two protons of the -CH2- linker are diastereotopic
    # and thus chemically non-equivalent.
    peaks_linker_ch2 = 2
    
    # b) The indazole moiety. This is a rigid, chiral structure with no internal symmetry plane.
    # We count the unique proton environments within it.
    
    # Proton on the pyrazole part of the indazole ring.
    peaks_indazole_pyrazole_h = 1
    
    # Protons of the three methyl groups on the camphor-derived skeleton. All are non-equivalent.
    peaks_indazole_methyls = 3
    
    # The single bridgehead methine (CH) proton.
    peaks_indazole_bridgehead_ch = 1
    
    # The four protons on the two aliphatic -CH2- groups of the bicyclic system. All are non-equivalent.
    peaks_indazole_aliphatic_ch2s = 4

    # Summing all the unique proton environments to get the total number of signals.
    total_peaks = (peaks_central_methyls + 
                   peaks_linker_ch2 + 
                   peaks_indazole_pyrazole_h + 
                   peaks_indazole_methyls + 
                   peaks_indazole_bridgehead_ch + 
                   peaks_indazole_aliphatic_ch2s)
                   
    # Print the explanation and final calculation.
    print("Based on the molecular structure, the expected number of 1H NMR peaks is calculated as follows:\n")
    print("1.  Central Core: The three methyl groups on the benzene ring are equivalent due to C3 symmetry.")
    print("    => 1 peak\n")
    print("2.  Substituent Arm (analyzing one of the three equivalent arms):")
    print("    a. Linker: The two protons on the -CH2- group are diastereotopic because the attached group is chiral.")
    print("       => 2 peaks")
    print("    b. Indazole Moiety: This chiral group has no internal symmetry.")
    print("       - One proton on the pyrazole ring => 1 peak")
    print("       - Three non-equivalent methyl groups => 3 peaks")
    print("       - One bridgehead methine (CH) proton => 1 peak")
    print("       - Four non-equivalent protons on the two aliphatic -CH2- groups => 4 peaks\n")
    
    print("Final Equation (sum of peaks from each non-equivalent set of protons):")
    # Using 'file=sys.stdout' to be explicit, which is good practice.
    print(f"{peaks_central_methyls} + {peaks_linker_ch2} + {peaks_indazole_pyrazole_h} + {peaks_indazole_methyls} + {peaks_indazole_bridgehead_ch} + {peaks_indazole_aliphatic_ch2s} = {total_peaks}", file=sys.stdout)


if __name__ == '__main__':
    calculate_nmr_peaks()
