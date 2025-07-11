import pandas as pd
from io import StringIO

def generate_xrd_indexing_table():
    """
    Identifies prominent peaks from the CuO XRD pattern and presents them
    in a table with their corresponding (h k l) Miller indices.
    """
    
    # Data derived from analyzing the provided XRD plot against standard CuO (tenorite) reference patterns (e.g., JCPDS No. 48-1548).
    # The 2-theta values are estimated from the plot.
    data = """
    Peak Position 2θ (degree), (h k l)
    ~35.5, (-1 1 1)
    ~38.7, (1 1 1)
    ~48.7, (-2 0 2)
    ~53.5, (0 2 0)
    ~58.3, (2 0 2)
    """
    
    # Using pandas for a well-formatted table output.
    # We read the data from the multiline string.
    s = StringIO(data)
    df = pd.read_csv(s, skipinitialspace=True)

    print("Indexing of Prominent XRD Peaks for CuO Nanoparticles")
    print("-----------------------------------------------------")
    # Using to_string() to make it look like a printed table
    print(df.to_string(index=False))

# Execute the function to print the table.
generate_xrd_indexing_table()
