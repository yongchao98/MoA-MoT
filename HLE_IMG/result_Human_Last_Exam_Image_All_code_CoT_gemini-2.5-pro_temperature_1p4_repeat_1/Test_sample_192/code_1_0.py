import pandas as pd

def index_cuo_xrd_peaks():
    """
    Identifies and indexes the prominent peaks from the provided CuO XRD spectrum.
    
    The peak positions are estimated from the provided graph for the pure CuO sample.
    The (h k l) Miller indices are assigned based on the standard diffraction data for
    monoclinic CuO (Tenorite), JCPDS card no. 48-1548.
    """
    
    # Data for the most prominent peaks estimated from the XRD plot
    # and indexed using standard CuO (Tenorite) data.
    # Format: (2-theta angle, Miller indices (h k l))
    peak_data = [
        (32.5, "(-1 1 0)"),
        (35.5, "(-1 1 1)"),
        (38.7, "(1 1 1)"),
        (48.7, "(-2 0 2)"),
        (53.5, "(0 2 0)"),
        (58.3, "(2 0 2)"),
        (61.5, "(-1 1 3)")
    ]

    # Create a DataFrame for better formatting
    df = pd.DataFrame(peak_data, columns=['Peak Position 2θ (degrees)', 'Miller Indices (h k l)'])

    # Print the table
    print("XRD Peak Indexing for CuO Nanoparticles")
    print("-" * 50)
    print(df.to_string(index=False))
    print("-" * 50)

# Run the function to generate and print the table
index_cuo_xrd_peaks()