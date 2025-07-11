import pandas as pd
try:
    from celegans_connectome.connectome import Connectome
except ImportError:
    print("Please install the required library by running: pip install c-elegans-connectome")
    exit()

def find_main_target_of_asj():
    """
    Finds the main projection target of ASJ axons based on synapse number
    in the C. elegans adult connectome.
    """
    try:
        # 1. Load the connectome data. This may download data on the first run.
        print("Loading C. elegans connectome data...")
        c = Connectome()
        connectome_df = c.get_connectome()
        print("Data loaded successfully.")

        # 2. Define the source neurons (ASJ left and right).
        source_neurons = ['ASJL', 'ASJR']

        # 3. Filter connections originating from ASJ neurons.
        asj_connections = connectome_df[connectome_df['pre'].isin(source_neurons)]

        # 4. Group by target cell and sum the number of synapses ('N').
        target_synapse_counts = asj_connections.groupby('post')['N'].sum()

        if target_synapse_counts.empty:
            print("No connections found for ASJ neurons.")
            return

        # 5. Find the target cell with the maximum number of synapses.
        main_target_cell = target_synapse_counts.idxmax()
        max_synapses = target_synapse_counts.max()

        # Find the individual contributions from ASJL and ASJR to the main target.
        synapses_from_asjl = asj_connections[
            (asj_connections['pre'] == 'ASJL') & (asj_connections['post'] == main_target_cell)
        ]['N'].sum()

        synapses_from_asjr = asj_connections[
            (asj_connections['pre'] == 'ASJR') & (asj_connections['post'] == main_target_cell)
        ]['N'].sum()

        # Print the final result and the equation.
        print("\n--- Analysis Result ---")
        print(f"The main projection target of ASJ axons is the cell: {main_target_cell}")
        print("\nThis is based on the total number of synapses received from the ASJ neuron pair (ASJL and ASJR).")
        print("The calculation is as follows:")
        print(f"Total Synapses = (Synapses from ASJL) + (Synapses from ASJR)")
        print(f"{max_synapses} = {synapses_from_asjl} + {synapses_from_asjr}")

    except Exception as e:
        print(f"\nAn error occurred: {e}")
        print("Please ensure you have a working internet connection for the initial data download.")
        print("You may need to install the library first: pip install c-elegans-connectome")

if __name__ == '__main__':
    find_main_target_of_asj()
