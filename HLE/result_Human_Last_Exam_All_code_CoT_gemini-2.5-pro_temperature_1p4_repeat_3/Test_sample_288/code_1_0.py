import pandas as pd
import requests
import io

def find_main_synaptic_target():
    """
    Finds the main synaptic target of a given source neuron in the C. elegans connectome.
    """
    # The URL points to the synaptic connectivity data from Cook et al. (2019)
    url = "https://www.cshlpress.com/supplemental/wormbook_218/wormbook_218_Cook_etal_supplemental_file_S2.csv"
    source_neurons = ['ASJL', 'ASJR']
    
    print(f"Analyzing connectome data to find the main projection target of {source_neurons}...\n")
    
    try:
        # Fetch the data from the URL
        response = requests.get(url)
        response.raise_for_status()  # This will raise an HTTPError for bad responses (4xx or 5xx)

        # Read the CSV data from the response content into a pandas DataFrame
        # Columns are: neuron_1 (presynaptic), neuron_2 (postsynaptic), type, n (number of synapses)
        connectome_df = pd.read_csv(io.StringIO(response.text))

        # Filter for rows where the presynaptic neuron is in our source_neurons list
        source_df = connectome_df[connectome_df['neuron_1'].isin(source_neurons)]

        # Filter for chemical synapses, which are denoted by 'S' (send)
        synapses_df = source_df[source_df['type'] == 'S'].copy()

        # Group by the target neuron ('neuron_2') and sum the number of synapses ('n')
        target_counts = synapses_df.groupby('neuron_2')['n'].sum()

        if target_counts.empty:
            print(f"No synaptic connections found from {source_neurons}.")
            return

        # Find the target cell with the maximum number of synapses
        main_target_cell = target_counts.idxmax()
        max_synapses = target_counts.max()

        # Print the final equation/result
        print(f"The number of synapses from ASJ to {main_target_cell} is {max_synapses}.")
        print(f"This is the highest number of synapses from ASJ to any single cell.")
        print(f"\nFinal Answer: The main projection target of ASJ axons is {main_target_cell}.")


    except requests.exceptions.RequestException as e:
        print(f"Error fetching data: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    find_main_synaptic_target()
