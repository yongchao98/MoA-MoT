import pandas as pd

def find_main_target():
    """
    Analyzes the C. elegans connectome to find the main synaptic target of ASJ neurons.
    """
    try:
        # The URL points to the N2U adult hermaphrodite connectome data (Witvliet et al., 2021)
        url = 'https://raw.githubusercontent.com/celegans/Connectome/master/Data/N2U/N2U_neuron_connections.csv'
        
        # Load the connectome data into a pandas DataFrame
        df = pd.read_csv(url)
        
        # ASJ is a neuron pair: ASJL and ASJR. Filter for outgoing connections from them.
        presynaptic_neurons = ['ASJL', 'ASJR']
        asj_connections = df[df['pre'].isin(presynaptic_neurons)]
        
        # We are interested in chemical synapses, counted in the 'N2U_S' column.
        # Group by the postsynaptic neuron ('post') and sum the synapse counts.
        synapse_counts_per_target = asj_connections.groupby('post')['N2U_S'].sum()
        
        if synapse_counts_per_target.empty:
            print("No synaptic targets found for ASJ neurons.")
            return

        # Find the target with the maximum number of synapses
        main_target_cell = synapse_counts_per_target.idxmax()
        
        # Get the individual connections from ASJL and ASJR to the main target
        main_target_connections = asj_connections[asj_connections['post'] == main_target_cell]
        
        print(f"The main projection target of ASJ axons is {main_target_cell}.")
        print("The calculation for the total number of synapses is as follows:")

        total_synapses = 0
        equation_parts = []
        
        # Iterate through the connections to the main target to build the output string
        for _, row in main_target_connections.iterrows():
            presynaptic_cell = row['pre']
            synapse_count = row['N2U_S']
            if synapse_count > 0:
                print(f"- Synapses from {presynaptic_cell} to {main_target_cell}: {synapse_count}")
                total_synapses += synapse_count
                equation_parts.append(str(synapse_count))

        equation_str = " + ".join(equation_parts)
        print(f"\nFinal Equation: {equation_str} = {total_synapses}")
        print(f"The cell with the highest synapse number from ASJ is {main_target_cell} with a total of {total_synapses} synapses.")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please check your internet connection and the data source URL.")

if __name__ == '__main__':
    find_main_target()