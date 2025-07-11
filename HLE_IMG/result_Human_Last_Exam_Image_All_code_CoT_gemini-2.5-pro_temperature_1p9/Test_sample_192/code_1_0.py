import pandas as pd
import io

def index_xrd_peaks():
    """
    This function identifies prominent peaks from the CuO XRD plot,
    assigns the standard Miller indices (hkl), and prints the results in a table.
    """
    # Data extracted from the XRD plot and correlated with standard CuO data (JCPDS 45-0937)
    # The peak at ~35.5 degrees is a combination of (-111) and (002) reflections,
    # and the peak at ~38.7 degrees combines (111) and (200) reflections.
    # We list the most commonly assigned indices.
    
    data = {
        "Observed Peak Position 2θ (degrees)": [32.5, 35.6, 38.8, 48.8, 53.5, 58.4],
        "Miller Indices (h k l)": ["(1 1 0)", "(-1 1 1)", "(1 1 1)", "(-2 0 2)", "(0 2 0)", "(2 0 2)"]
    }

    # Create a DataFrame using pandas
    df = pd.DataFrame(data)

    # Print the DataFrame as a nicely formatted table
    print("XRD Peak Indexing for CuO Nanoparticles")
    print("-" * 50)
    print(df.to_string(index=False))
    print("-" * 50)
    print("Note: Peak positions are estimated from the provided graph.")
    print("Miller indices (h k l) are assigned based on the standard pattern for monoclinic CuO (Tenorite).")

# Execute the function to display the table
index_xrd_peaks()