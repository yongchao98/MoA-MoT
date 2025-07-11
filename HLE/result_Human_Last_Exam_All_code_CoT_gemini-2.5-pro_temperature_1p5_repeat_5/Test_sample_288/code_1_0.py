import pandas as pd
import requests
import io

def find_asj_main_target():
    """
    This script identifies the main projection target of the ASJ neurons in the
    adult C. elegans hermaphrodite connectome by analyzing public data.
    """
    try:
        # The URL points to the supplementary data from Cook et al., Nature, 2019
        url = "http://wormwiring.org/si/SI-1-Connectome.xlsx"
        print(f"Fetching C. elegans connectome data from {url}...")
        # Use a timeout for the request to prevent indefinite hanging
        response = requests.get(url, timeout=30)
        # Raise an HTTPError for bad responses (4xx or 5xx)
        response.raise_for_status()
        print("Data downloaded successfully.\n")

        # Read the Excel file content from memory
        excel_file = io.BytesIO(response.content)
        df = pd.read_excel(excel_file, sheet_name='Connectome')

    except requests.exceptions.RequestException as e:
        print(f"Error fetching data: {e}")
        print("Failed to retrieve connectome data. Please check your internet connection or if the data URL is still valid.")
        return
    except Exception as e:
        print(f"An error occurred while processing the data file: {e}")
        return

    # Filter for connections where the presynaptic neuron is ASJL or ASJR
    asj_connections = df[df['Neuron 1'].isin(['ASJL', 'ASJR'])].copy()

    def get_neuron_class(neuron_name):
        """Removes the trailing 'L' or 'R' to identify the neuron class."""
        if isinstance(neuron_name, str) and (neuron_name.endswith('L') or neuron_name.endswith('R')):
            return neuron_name[:-1]
        return neuron_name

    # Create a new column for the target neuron's class
    asj_connections['Target Class'] = asj_connections['Neuron 2'].apply(get_neuron_class)

    # Group by the target class and sum the number of synapses ('N')
    synapse_sums = asj_connections.groupby('Target Class')['N'].sum()

    if synapse_sums.empty:
        print("No synaptic connections were found for ASJ neurons in the dataset.")
        return

    # Identify the class with the maximum number of synapses
    main_target_class = synapse_sums.idxmax()
    total_synapses = int(synapse_sums.max())

    print(f"The main projection target of ASJ axons is the {main_target_class} neuron class.")
    print("--------------------------------------------------")
    print(f"Breakdown of synapses from ASJ to the {main_target_class} class:")

    # Filter connections to the main target class to show the detailed breakdown
    target_connections = asj_connections[asj_connections['Target Class'] == main_target_class]
    
    equation_parts = []
    # Sort for consistent output
    for _, row in target_connections.sort_values(by=['Neuron 1', 'Neuron 2']).iterrows():
        pre_neuron = row['Neuron 1']
        post_neuron = row['Neuron 2']
        synapse_count = row['N']
        print(f"{pre_neuron} -> {post_neuron}: {synapse_count} synapses")
        equation_parts.append(str(synapse_count))

    print("--------------------------------------------------")
    equation_str = " + ".join(equation_parts)
    print(f"Total synapses = {equation_str} = {total_synapses}")
    
    # The final answer in the required format
    print(f"\n<<< {main_target_class} >>>")

if __name__ == "__main__":
    find_asj_main_target()