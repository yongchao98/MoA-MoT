import pandas as pd

def find_main_asj_target():
    """
    Finds the main projection target of ASJ axons in the C. elegans connectome
    based on the number of synapses.
    """
    try:
        # URL to the N2U adult C. elegans connectome data
        url = "https://raw.githubusercontent.com/Vassar-CIL/c-elegans-connectome/master/data/N2U_neuron_connections.csv"
        
        # Read the connectome data into a pandas DataFrame
        connectome_df = pd.read_csv(url)
        
        # Filter for connections where ASJ is the presynaptic neuron
        asj_projections = connectome_df[connectome_df['pre'] == 'ASJ']
        
        # Group by the postsynaptic neuron ('post') and sum the number of synapses ('N')
        # This accounts for multiple types of connections (e.g., chemical, electrical) to the same target
        target_synapse_counts = asj_projections.groupby('post')['N'].sum()
        
        # Find the target cell with the maximum number of synapses
        main_target_cell = target_synapse_counts.idxmax()
        max_synapses = target_synapse_counts.max()
        
        print(f"The main projection target of ASJ axons is the cell '{main_target_cell}'.")
        print(f"Total number of synapses from ASJ to {main_target_cell}: {max_synapses}")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Could not fetch or process the connectome data.")

if __name__ == "__main__":
    find_main_asj_target()