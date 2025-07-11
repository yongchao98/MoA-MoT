import pandas as pd
from io import StringIO

def generate_xrd_index_table():
    """
    Identifies and indexes the prominent peaks from an XRD pattern of CuO nanoparticles
    and prints the results in a formatted table.
    The indexing is based on the standard data for monoclinic CuO (Tenorite),
    JCPDS card no. 45-0937.
    """
    
    # Data derived from the graph and standard CuO XRD patterns
    # Columns: Peak Position (2θ), Miller Indices (hkl)
    # The two most prominent peaks are at ~35.5 and ~38.7 degrees.
    peak_data = {
        'Peak Position 2θ (degree)': ["~32.5", "~35.5", "~38.7", "~48.7", "~53.5", "~58.3"],
        'Miller Indices (h k l)': ["(1 1 0)", "(-1 1 1)", "(1 1 1)", "(-2 0 2)", "(0 2 0)", "(2 0 2)"]
    }

    # Create a pandas DataFrame for nice formatting
    df = pd.DataFrame(peak_data)

    print("Table of Indexed XRD Peaks for CuO Nanoparticles")
    print("-" * 52)
    print(df.to_string(index=False))
    print("-" * 52)
    print("\nNote: Peak positions are estimated from the provided graph and are consistent")
    print("with the standard pattern for monoclinic CuO (Tenorite).")

# Execute the function to print the table
generate_xrd_index_table()