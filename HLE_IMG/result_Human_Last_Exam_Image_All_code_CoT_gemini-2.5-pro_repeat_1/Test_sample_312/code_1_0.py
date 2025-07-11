import collections

def solve():
    """
    This function solves the SWNT plot identification problem and prints the result.
    
    The logic for identifying the m-value for each plot is as follows:
    
    1.  Properties of (4, m) SWNTs:
        - (4,0): Semiconducting, Zigzag (achiral), Chiral Angle=0°, Smallest diameter among semiconductors.
        - (4,1): Metallic, Chiral, Chiral Angle=10.9°.
        - (4,2): Semiconducting, Chiral, Chiral Angle=19.1°, Medium diameter.
        - (4,3): Semiconducting, Chiral, Chiral Angle=25.3°, Largest diameter.
        - (4,4): Metallic, Armchair (achiral), Chiral Angle=30°.

    2.  Analysis of Oscillator Strength Plots (#2, #4, #6, #8):
        - The tilt angle of the cutting lines increases with the chiral angle.
        - Plot #2: Vertical lines (0° angle) -> m=0.
        - Plot #8: Smallest tilt -> m=1.
        - Plot #4: Medium tilt -> m=2.
        - Plot #6: Largest tilt -> m=3.
        - The plot for m=4 (Armchair) is the one that is omitted.

    3.  Analysis of Dipole Moment (Band Structure) Plots (#1, #3, #5, #7, #9):
        - Metallic (no band gap): Plots #1 and #5.
            - Plot #1: High symmetry, fewer bands. Characteristic of Armchair (4,4) -> m=4.
            - Plot #5: Less symmetric, more bands. Must be the other metallic tube, Chiral (4,1) -> m=1.
        - Semiconducting (has a band gap): Plots #3, #7, #9.
            - The band gap is inversely proportional to the diameter. Diameter order: (4,0) < (4,2) < (4,3).
            - So, band gap order: (4,0) > (4,2) > (4,3).
            - The number of bands increases with complexity/diameter.
            - Plot #9: Largest gap, fewest bands -> m=0.
            - Plot #7: Medium gap, medium bands -> m=2.
            - Plot #3: Smallest gap, most bands -> m=3.

    4.  Final Sequence Assembly:
        Combining the results for each plot number from 1 to 9.
    """
    
    # m-value for each plot from 1 to 9
    # Plot #1: m=4
    # Plot #2: m=0
    # Plot #3: m=3
    # Plot #4: m=2
    # Plot #5: m=1
    # Plot #6: m=3
    # Plot #7: m=2
    # Plot #8: m=1
    # Plot #9: m=0
    m_values = [4, 0, 3, 2, 1, 3, 2, 1, 0]
    
    # Format the output string as requested
    output_string = "{" + ", ".join(map(str, m_values)) + "}"
    
    print(output_string)

solve()