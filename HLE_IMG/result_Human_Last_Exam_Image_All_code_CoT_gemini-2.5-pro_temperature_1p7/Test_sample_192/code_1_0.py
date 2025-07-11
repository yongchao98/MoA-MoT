import pandas as pd

def index_xrd_peaks():
    """
    Identifies and indexes the prominent XRD peaks for CuO nanoparticles
    and prints them in a formatted table.
    
    The peak positions are estimated from the provided graph, and the Miller indices (hkl)
    are assigned based on the standard diffraction pattern of monoclinic CuO (JCPDS card no. 45-0937).
    """
    
    # Data derived from visual inspection of the XRD plot and comparison with standard CuO data
    peak_data = {
        'Peak Position 2θ (degree)': [
            "~32.5",
            "~35.5",
            "~38.7",
            "~48.7",
            "~53.5",
            "~58.3",
            "~61.5"
        ],
        'Miller Indices (h k l)': [
            "(110)",
            "(002), (-111)",
            "(111), (200)",
            "(-202)",
            "(020)",
            "(202)",
            "(-113)"
        ]
    }
    
    # Create a pandas DataFrame for better formatting
    df = pd.DataFrame(peak_data)
    
    # Print the table
    print("XRD Peak Indexing for CuO Nanoparticles")
    print("-" * 50)
    print(df.to_string(index=False))
    print("-" * 50)

# Execute the function to generate the table
index_xrd_peaks()