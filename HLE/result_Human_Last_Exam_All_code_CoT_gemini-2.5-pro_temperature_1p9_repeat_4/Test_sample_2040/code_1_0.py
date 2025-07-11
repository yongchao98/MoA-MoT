#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def analyze_spdc_in_borophene():
    """
    Analyzes whether free-standing boron nanosheets (borophene) would
    exhibit Spontaneous Parametric Down-Conversion (SPDC) and prints the
    reasoning.
    """

    print("Question: Would free-standing boron nanosheets be expected to exhibit spontaneous parametric downconversion (SPDC)?")
    print("-" * 80)

    # Step 1: Explain the physical requirement for SPDC
    print("\nStep 1: The Physics of Spontaneous Parametric Down-Conversion (SPDC)")
    print("  - SPDC is a 'second-order' nonlinear optical process where one photon splits into two.")
    print("  - A material's ability to facilitate this is measured by its second-order nonlinear susceptibility, denoted as χ⁽²⁾.")
    print("  - A core requirement for a material to exhibit SPDC is to have a non-zero value for χ⁽²⁾.")
    
    # Step 2: Analyze the relevant properties of boron nanosheets
    print("\nStep 2: Material Properties of Boron Nanosheets (Borophene)")
    print("  - A material's χ⁽²⁾ value is determined by its crystal symmetry.")
    print("  - In materials that have a center of inversion symmetry (i.e., they are 'centrosymmetric'), the bulk second-order susceptibility is zero.")
    print("  - The equation for such materials is:")
    chi_2 = 0
    print(f"      χ⁽²⁾ = {chi_2}")
    print("  - The common and most stable phases of free-standing borophene are known to be centrosymmetric.")

    # Step 3: Combine the facts to form a conclusion
    print("\nStep 3: Conclusion")
    print("  - Since borophene is generally a centrosymmetric material, its bulk second-order susceptibility (χ⁽²⁾) is expected to be zero.")
    print("  - Without a non-zero χ⁽²⁾, the material cannot support the SPDC process.")
    print("  - Therefore, free-standing boron nanosheets are not expected to exhibit SPDC.")
    print("\n  (Note: Symmetry can be broken at material edges or by defects, which could lead to a very weak, localized effect, but it is not an intrinsic property of the bulk material.)")


if __name__ == "__main__":
    analyze_spdc_in_borophene()
