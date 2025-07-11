import pandas as pd

def solve_microscopy_problem():
    """
    Solves the problem by identifying all fluorophores, their excitation
    wavelengths, and determining the first one to be used based on imaging best practices.
    """

    # Step 1 & 2: Identify fluorophores and their excitation wavelengths.
    # The transgenic zebrafish has eGFP and DsRed.
    # A far-red HaloTag ligand (similar to JF646) is added.
    fluorophores = {
        'eGFP in Neutrophils': 488,
        'DsRed in Macrophages': 559,
        'Far-Red Dye on HaloTag in Neutrophils': 630
    }

    print("Step 1: Identify all fluorescent sources in the sample.")
    print("The sample contains a transgenic zebrafish with three potential fluorescent signals:")
    for component, wavelength in fluorophores.items():
        print(f"- {component} (Excitation: ~{wavelength} nm)")
    print("-" * 30)

    # Step 3: Interpret "first" based on imaging best practices.
    # Best practice for live, multi-color imaging is to start with the longest
    # wavelength (lowest energy) to minimize phototoxicity and photobleaching of other dyes.
    wavelengths = list(fluorophores.values())
    
    # Sort wavelengths in descending order to find the first one to use.
    wavelengths.sort(reverse=True)
    
    first_wavelength = wavelengths[0]

    print("Step 2: Determine the order of imaging.")
    print("In multi-color live imaging, the standard procedure is to excite fluorophores sequentially from the longest to the shortest wavelength.")
    print("This minimizes phototoxicity and prevents photobleaching of sensitive dyes by high-energy light.")
    print("\nComparing the excitation wavelengths:")
    
    # "Output each number in the final equation!" -> printing the comparison
    comparison_str = ' > '.join(map(str, wavelengths))
    print(f"The imaging order should be based on the principle: {comparison_str}")

    print("\nStep 3: Conclusion.")
    print(f"Therefore, the first excitation wavelength that should be used in the experiment is the longest one, which is {first_wavelength} nm.")

solve_microscopy_problem()