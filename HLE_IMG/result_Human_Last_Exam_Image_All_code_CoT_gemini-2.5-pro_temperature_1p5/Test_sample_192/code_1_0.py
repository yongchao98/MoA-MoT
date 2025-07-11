import pandas as pd

def generate_peak_table():
    """
    This function generates and prints a table of XRD peak positions and their
    corresponding Miller indices (h k l) for the provided CuO nanoparticle data.
    """
    # Data derived from analyzing the XRD plot against the standard JCPDS card for CuO (e.g., #05-0661).
    # Peak positions are visually estimated from the provided graph.
    peak_data = {
        'Peak Position (2θ°)': ['~32.5', '~35.5', '~38.7', '~48.7', '~53.5', '~58.3'],
        'Miller Indices (h k l)': ['(1 1 0)', '(-1 1 1)', '(1 1 1)', '(-2 0 2)', '(0 2 0)', '(2 0 2)']
    }

    # Create a pandas DataFrame for nice formatting.
    df = pd.DataFrame(peak_data)

    # Print the table.
    print("XRD Peak Indexing for CuO Nanoparticles")
    print("-" * 45)
    print(df.to_string(index=False))
    print("-" * 45)
    print("\nNote: Peak positions are approximate values estimated from the graph.")

# Run the function to display the table.
generate_peak_table()