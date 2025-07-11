import pandas as pd

def index_xrd_peaks():
    """
    Identifies prominent XRD peaks for CuO from the provided plot,
    indexes them with Miller indices (h k l), and prints the result in a table.
    The indexing is based on the standard data for monoclinic CuO (Tenorite).
    """

    # Data extracted from the plot and standard CuO diffraction patterns (JCPDS No. 45-0937)
    # The two most prominent peaks are at ~35.5 and ~38.7 degrees.
    peak_data = {
        'Peak Position 2θ (degree)': [32.5, 35.5, 38.7, 48.7, 53.5, 58.3],
        'Miller Index (h k l)': ['(1 1 0)', '(-1 1 1)', '(1 1 1)', '(-2 0 2)', '(0 2 0)', '(2 0 2)']
    }

    # Create a pandas DataFrame for a well-formatted table
    df = pd.DataFrame(peak_data)

    print("XRD Peak Indexing for CuO Nanoparticles")
    print("-" * 50)
    # Print the DataFrame without the pandas index
    print(df.to_string(index=False))
    print("-" * 50)
    print("\nNote:")
    print("1. Peak positions (2θ) are estimated from the provided graph.")
    print("2. Miller indices (h k l) correspond to the monoclinic structure of CuO.")


# Execute the function
index_xrd_peaks()