import pandas as pd

def index_xrd_peaks():
    """
    Identifies prominent peaks from the CuO XRD pattern and presents them with their
    corresponding Miller indices (hkl) in a table.
    The indexing is based on the standard monoclinic structure of CuO (tenorite).
    """

    # Data for prominent peaks identified from the XRD plot
    # Peak positions (2-theta) are estimated from the provided graph.
    # (hkl) values are from standard CuO diffraction data (JCPDS card no. 45-0937).
    peak_data = {
        'Peak Position 2θ (degree)': [
            "~32.5",
            "~35.5",
            "~38.7",
            "~48.7",
            "~53.4",
            "~58.3"
        ],
        'Miller Indices (h k l)': [
            "(1 1 0)",
            "(0 0 2) / (1 1 -1)",
            "(1 1 1)",
            "(-2 0 2)",
            "(0 2 0)",
            "(2 0 2)"
        ]
    }

    # Create a pandas DataFrame for neat formatting
    df = pd.DataFrame(peak_data)

    # Print the table
    print("XRD Peak Indexing for Monoclinic CuO Nanoparticles")
    print("-" * 55)
    print(df.to_string(index=False))
    print("-" * 55)
    print("\nNote: Peak positions are approximate values estimated from the provided graph.")

if __name__ == "__main__":
    index_xrd_peaks()