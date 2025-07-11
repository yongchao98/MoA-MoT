import pandas as pd

def generate_xrd_index_table():
    """
    Identifies prominent peaks from the CuO XRD plot, assigns Miller indices (hkl),
    and prints the results in a formatted table.
    """

    # Data for the prominent peaks identified from the XRD pattern.
    # Peak positions (2-theta in degrees) are estimated from the graph.
    # Miller indices (hkl) are assigned based on the standard JCPDS card no. 45-0937 for monoclinic CuO.
    data = {
        'Observed 2θ (degree)': [
            '~32.5',
            '~35.5',
            '~38.7',
            '~48.7',
            '~53.5',
            '~58.3'
        ],
        '(h k l)': [
            '(1 1 0)',
            '(-1 1 1)',
            '(1 1 1)',
            '(-2 0 2)',
            '(0 2 0)',
            '(2 0 2)'
        ]
    }

    # Create a DataFrame for better formatting
    df = pd.DataFrame(data)
    df.index = df.index + 1  # Start peak numbering from 1
    df.index.name = "Peak No."

    # Print the header and the formatted table
    print("Indexing of Prominent XRD Peaks for CuO Nanoparticles")
    print("-" * 55)
    print(df.to_string())
    print("-" * 55)


if __name__ == '__main__':
    generate_xrd_index_table()
