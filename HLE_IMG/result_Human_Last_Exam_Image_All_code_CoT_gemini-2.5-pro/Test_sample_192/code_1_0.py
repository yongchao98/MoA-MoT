import pandas as pd
from io import StringIO

def generate_xrd_table():
    """
    Generates and prints a table of indexed XRD peaks for CuO nanoparticles.
    The peak positions are estimated from the provided graph and indexed
    based on the standard pattern for monoclinic CuO (Tenorite).
    """

    # Data for the prominent peaks
    data = """2θ (approx. degree),(h k l) Index
32.5,(110)
35.5,(-111)
38.7,(111)
48.7,(-202)
53.5,(020)
58.3,(202)
"""
    
    # Read the data into a pandas DataFrame
    df = pd.read_csv(StringIO(data))

    # Print the table
    print("Indexed XRD Peaks for CuO Nanoparticles")
    print("="*40)
    print(df.to_string(index=False))
    print("="*40)
    print("\nNote: Peak positions are estimated from the graph.")
    print("Indexing corresponds to standard monoclinic CuO.")

if __name__ == '__main__':
    generate_xrd_table()
