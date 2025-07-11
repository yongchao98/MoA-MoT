import collections

def solve_nanotube_puzzle():
    """
    This script determines the chiral index 'm' for each of the nine plots
    based on the physical properties of single-walled carbon nanotubes (SWNTs).
    """

    # Step 1: Define chiralities and their properties
    # n is fixed at 4. m varies from 0 to 4.
    # A nanotube is metallic if (n-m) is a multiple of 3.
    # Diameter is proportional to sqrt(n^2 + n*m + m^2). Band gap is inversely proportional to diameter.
    chiralities = {
        0: {'type': 'Semiconducting (Zigzag)', 'n-m': 4, 'gap_rank': 1, 'plot_pair': (9, 6)},
        1: {'type': 'Metallic (Chiral)', 'n-m': 3, 'gap_rank': None, 'plot_pair': (5, 4)},
        2: {'type': 'Semiconducting (Chiral)', 'n-m': 2, 'gap_rank': 2, 'plot_pair': (7, 8)},
        3: {'type': 'Semiconducting (Chiral)', 'n-m': 1, 'gap_rank': 3, 'plot_pair': (3, 2)},
        4: {'type': 'Metallic (Armchair)', 'n-m': 0, 'gap_rank': None, 'plot_pair': (1, None)} # Oscillator plot is omitted
    }

    # Step 2: Map plots to m-values based on analysis
    # The analysis pairs each (n,m) value with a band structure plot and an oscillator strength plot.
    # For example, (4,0) has the largest band gap (Plot #9) and vertical cutting lines (Plot #6).
    plot_to_m_mapping = collections.OrderedDict()
    for m_val, properties in chiralities.items():
        band_plot, oscillator_plot = properties['plot_pair']
        if band_plot is not None:
            plot_to_m_mapping[band_plot] = m_val
        if oscillator_plot is not None:
            plot_to_m_mapping[oscillator_plot] = m_val
    
    # Sort by plot number to get the final sequence
    sorted_mapping = collections.OrderedDict(sorted(plot_to_m_mapping.items()))

    print("Assignment of m-value to each plot:")
    final_sequence = []
    for plot_num, m_val in sorted_mapping.items():
        print(f"Plot #{plot_num} corresponds to m = {m_val}")
        final_sequence.append(m_val)

    # Step 3: Format and print the final result
    # The final answer is the sequence of m-values for plots 1 through 9.
    final_sequence_str = ", ".join(map(str, final_sequence))
    
    print("\nThe final sequence of m-values for plots #1 through #9 is:")
    print(f"{{{final_sequence_str}}}")

solve_nanotube_puzzle()