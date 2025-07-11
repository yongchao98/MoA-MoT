import pandas as pd

def index_xrd_peaks():
    """
    Identifies and indexes the prominent XRD peaks for CuO nanoparticles
    based on the provided plot and standard crystallographic data for
    monoclinic CuO (tenorite). It then prints the results in a table.
    """
    
    # Data extracted from the XRD plot and correlated with standard CuO patterns.
    # The data contains the approximate 2-theta angle and the corresponding
    # Miller indices (h k l) for the crystal planes.
    data = {
        'Observed Peak Position (2θ, approx.)': [
            '32.5°', 
            '35.5°', 
            '38.7°', 
            '48.8°', 
            '53.5°', 
            '58.3°'
        ],
        'Miller Indices (h k l)': [
            '(1 1 0)',
            '(-1 1 1) / (0 0 2)',
            '(1 1 1) / (2 0 0)',
            '(-2 0 2)',
            '(0 2 0)',
            '(2 0 2)'
        ]
    }
    
    # Create a pandas DataFrame for better formatting
    df = pd.DataFrame(data)
    
    # Print the DataFrame as a string, which gives a nicely formatted table
    print("--- XRD Peak Indexing for CuO Nanoparticles ---")
    print(df.to_string(index=False))
    print("\nNote: Peak positions are estimated from the provided plot. Indexing is")
    print("based on the standard pattern for monoclinic CuO (Tenorite). Some observed")
    print("peaks correspond to multiple, closely spaced diffraction planes.")

# Run the function to display the results
index_xrd_peaks()