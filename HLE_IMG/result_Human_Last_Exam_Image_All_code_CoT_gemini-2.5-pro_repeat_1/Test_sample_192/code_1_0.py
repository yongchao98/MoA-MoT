import pandas as pd

def index_cuo_xrd_peaks():
    """
    This function identifies prominent XRD peaks for CuO nanoparticles from the provided
    graph and assigns the corresponding Miller indices (h k l) based on standard data
    for the monoclinic CuO structure (Tenorite).

    The output is a formatted table of the results.
    """
    # Data is based on visual estimation from the graph and comparison with
    # standard JCPDS card no. 48-1548 for monoclinic CuO.
    peak_data = {
        'Approx. 2θ (degree)': [32.5, 35.5, 38.7, 48.8, 53.5, 58.3, 61.5],
        '(h k l) Index': ['(110)', '(-111)', '(111)', '(-202)', '(020)', '(202)', '(-113)']
    }

    # Create a pandas DataFrame for neat printing
    df = pd.DataFrame(peak_data)

    # Print the DataFrame as a string, which provides a clean table format
    print("XRD Peak Indexing for CuO Nanoparticles")
    print("-" * 40)
    print(df.to_string(index=False))
    print("-" * 40)
    print("Note: Peak positions are estimated from the provided XRD plot.")

# Execute the function to print the table
index_cuo_xrd_peaks()