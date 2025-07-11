import pandas as pd

def generate_xrd_table():
    """
    Generates and prints a table of indexed XRD peaks for CuO nanoparticles.
    """
    # Data is based on visual inspection of the plot and comparison with
    # standard CuO (Tenorite, e.g., JCPDS card no. 48-1548) XRD data.
    # Each entry is a dictionary representing a peak.
    peak_data = [
        {'Approximate 2θ (degree)': 32.5, 'Miller Indices (h k l)': '(1 1 0)'},
        {'Approximate 2θ (degree)': 35.5, 'Miller Indices (h k l)': '(-1 1 1)'},
        {'Approximate 2θ (degree)': 38.7, 'Miller Indices (h k l)': '(1 1 1)'},
        {'Approximate 2θ (degree)': 48.8, 'Miller Indices (h k l)': '(-2 0 2)'},
        {'Approximate 2θ (degree)': 53.5, 'Miller Indices (h k l)': '(0 2 0)'},
        {'Approximate 2θ (degree)': 58.3, 'Miller Indices (h k l)': '(2 0 2)'},
        {'Approximate 2θ (degree)': 61.5, 'Miller Indices (h k l)': '(-1 1 3)'}
    ]

    # Create a pandas DataFrame for better formatting
    df = pd.DataFrame(peak_data)

    # Print the table
    print("Indexing of Prominent XRD Peaks for CuO Nanoparticles")
    print(df.to_string(index=False))

if __name__ == '__main__':
    generate_xrd_table()