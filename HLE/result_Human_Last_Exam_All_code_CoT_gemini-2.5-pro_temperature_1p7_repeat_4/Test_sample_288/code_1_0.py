import pandas as pd
import io
import requests

def find_main_synaptic_target():
    """
    This script finds the main synaptic target of the ASJ neurons in the
    adult C. elegans hermaphrodite connectome.
    """
    # The URL points to the synapse list from the Witvliet et al. 2021 dataset
    # for the N2U adult hermaphrodite.
    url = 'http://wormwiring.org/si/N2U_hermaphrodite_synapses.csv'

    try:
        # Fetch the data from the URL
        response = requests.get(url)
        response.raise_for_status()  # This will raise an error for a bad response (404, 500, etc.)

        # Read the CSV data from the response text into a pandas DataFrame
        # The columns are: cell_pre, cell_post, N, type
        synapse_data = pd.read_csv(io.StringIO(response.text))

        # Define the source neurons of interest (ASJ left and right)
        source_neurons = ['ASJL', 'ASJR']

        # Filter the data for chemical synapses originating from ASJ neurons
        asj_outputs = synapse_data[
            (synapse_data['cell_pre'].isin(source_neurons)) &
            (synapse_data['type'] == 'chemical')
        ]

        # Group the filtered data by the postsynaptic cell ('cell_post') and
        # sum the number of synapses ('N') for each target cell.
        target_synapse_counts = asj_outputs.groupby('cell_post')['N'].sum()

        if target_synapse_counts.empty:
            print("No chemical synapse targets found for ASJ neurons.")
            return

        # Find the target cell with the maximum total number of synapses
        main_target_cell = target_synapse_counts.idxmax()
        max_synapses = target_synapse_counts.max()

        # To show the equation, we calculate the contribution from each ASJ neuron
        synapses_from_asjl = asj_outputs[
            (asj_outputs['cell_pre'] == 'ASJL') & (asj_outputs['cell_post'] == main_target_cell)
        ]['N'].sum()

        synapses_from_asjr = asj_outputs[
            (asj_outputs['cell_pre'] == 'ASJR') & (asj_outputs['cell_post'] == main_target_cell)
        ]['N'].sum()


        print(f"The main projection target of ASJ axons is {main_target_cell}, which receives a total of {int(max_synapses)} synapses.")
        print("The calculation is as follows:")
        # Print the final equation with each number
        print(f"{int(synapses_from_asjl)} (from ASJL) + {int(synapses_from_asjr)} (from ASJR) = {int(max_synapses)}")


    except requests.exceptions.RequestException as e:
        print(f"Failed to fetch data from URL. Error: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    find_main_synaptic_target()
