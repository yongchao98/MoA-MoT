import pandas as pd
import io
import requests

def find_main_synapse_target():
    """
    This function finds the main projection target of a specified neuron class
    in the C. elegans connectome based on the number of synapses.
    """
    # URL for the N2U adult hermaphrodite chemical synapse data (Cook et al., 2019)
    url = "https://docs.google.com/spreadsheets/d/e/2PACX-1vRk_91nJ7_BH_6A2ssI2x5X3_dCslxS25y3SCp7Sy2W1242n7i52Byl2-wH3zLp_VTr2D35l2jV8hX5/pub?gid=2040523384&single=true&output=csv"
    
    source_neurons = ['ASJL', 'ASJR']

    try:
        # Fetch the data from the URL
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for HTTP errors
        
        # Read the CSV data into a pandas DataFrame
        connectome_df = pd.read_csv(io.StringIO(response.text))
        connectome_df.columns = ['Source', 'Target', 'SynapseCount']
        
        # Filter for connections where the source is in our list of source neurons
        source_outputs = connectome_df[connectome_df['Source'].isin(source_neurons)]
        
        if source_outputs.empty:
            print(f"No connections found from neurons: {', '.join(source_neurons)}")
            return
            
        # Group by the target neuron and sum the synapse counts
        target_synapse_counts = source_outputs.groupby('Target')['SynapseCount'].sum()
        
        # Find the target with the maximum number of synapses
        main_target_cell = target_synapse_counts.idxmax()
        total_synapses = target_synapse_counts.max()
        
        # Get the individual synapse counts from each source neuron to the main target
        breakdown = source_outputs[source_outputs['Target'] == main_target_cell]
        
        print(f"The main projection target of ASJ axons is '{main_target_cell}'.")
        print("\nThis is calculated by summing the synapses from ASJL and ASJR to the target cell:")
        
        equation_parts = []
        for index, row in breakdown.iterrows():
            print(f"Synapses from {row['Source']} to {row['Target']}: {row['SynapseCount']}")
            equation_parts.append(str(row['SynapseCount']))
            
        final_equation = " + ".join(equation_parts)
        
        print(f"\nFinal equation: {final_equation} = {total_synapses}")

    except requests.exceptions.RequestException as e:
        print(f"Error: Could not retrieve data from the network. Details: {e}")
    except Exception as e:
        print(f"An error occurred during processing: {e}")

if __name__ == "__main__":
    find_main_synapse_target()