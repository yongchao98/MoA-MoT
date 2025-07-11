import pandas as pd

def index_xrd_peaks():
    """
    Identifies and indexes the prominent peaks from the provided CuO XRD pattern.
    The peak positions are estimated from the graph, and the (h k l) indices
    are based on the standard monoclinic structure of CuO (Tenorite).
    """

    # Data for prominent peaks: (2-theta value, (hkl) index)
    # The two most intense peaks are at ~35.6 and ~38.8 degrees.
    peak_data = {
        'Peak Position (2θ, degree)': [32.5, 35.6, 38.8, 48.8, 53.5, 58.4],
        'Miller Indices (h k l)': ["(110)", "(-111)", "(111)", "(-202)", "(020)", "(202)"]
    }

    # Create a pandas DataFrame for better formatting
    df = pd.DataFrame(peak_data)

    # Print the table
    print("Indexing of Prominent XRD Peaks for CuO Nanoparticles")
    print("="*55)
    print(df.to_string(index=False))
    print("="*55)
    print("\nNote:")
    print("1. Peak positions are estimated from the provided XRD plot for the pure CuO sample.")
    print("2. The Miller indices (h k l) correspond to the standard monoclinic crystal structure of CuO (Tenorite).")

if __name__ == '__main__':
    index_xrd_peaks()