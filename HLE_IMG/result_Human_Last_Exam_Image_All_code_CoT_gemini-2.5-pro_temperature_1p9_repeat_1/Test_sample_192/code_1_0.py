import pandas as pd
from io import StringIO

def generate_xrd_index_table():
    """
    Analyzes the prominent peaks from the provided CuO XRD plot and indexes them
    based on standard crystallographic data for monoclinic CuO.
    
    The function then prints a formatted table of the indexed peaks.
    """
    
    # Data derived from analyzing the plot and comparing with standard JCPDS data for monoclinic CuO.
    # The most prominent peaks are selected for indexing.
    # The peak positions are estimated from the provided image.
    data = """
    2θ (degree), (h k l)
    32.5, (1 1 0)
    35.6, (-1 1 1)
    38.7, (1 1 1)
    48.8, (-2 0 2)
    53.5, (0 2 0)
    58.3, (2 0 2)
    61.5, (-1 1 3)
    """
    
    # Read the data into a pandas DataFrame for nice formatting
    df = pd.read_csv(StringIO(data), sep=',', skipinitialspace=True)
    
    print("Indexing of Prominent XRD Peaks for CuO Nanoparticles")
    print("=" * 55)
    print(df.to_string(index=False))
    print("=" * 55)
    print("\nNote:")
    print("1. Peak positions (2θ) are estimated from the provided plot.")
    print("2. (h k l) indices are based on the standard reference pattern for")
    print("   monoclinic CuO (JCPDS card no. 45-0937).")

if __name__ == '__main__':
    generate_xrd_index_table()