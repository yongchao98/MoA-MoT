import pandas as pd
import io
import requests

def find_asj_main_target():
    """
    Finds the main synaptic target of the ASJ neuron in C. elegans.

    This function fetches the C. elegans connectome data, filters for presynaptic
    connections from ASJ, and identifies the postsynaptic neuron with the
    highest number of synapses.
    """
    # URL for the adult C. elegans hermaphrodite connectome data
    url = "http://www.wormwiring.org/db/DATA/herm_full_edgelist.csv"

    print(f"Fetching connectome data from {url}...")

    try:
        # Get the data from the URL
        response = requests.get(url)
        response.raise_for_status()  # Check for HTTP errors

        # The CSV data is preceded by comments. We need to find the header row.
        content = response.text
        lines = content.splitlines()
        
        # Find the line index where the actual CSV header starts
        header_line_index = -1
        for i, line in enumerate(lines):
            if line.strip().startswith('neuron_1,neuron_2,type,N'):
                header_line_index = i
                break
        
        if header_line_index == -1:
            print("Error: Could not find the CSV header in the downloaded data.")
            return

        # Create a clean CSV string starting from the header
        csv_data = "\n".join(lines[header_line_index:])

        # Read the cleaned CSV data into a pandas DataFrame
        df = pd.read_csv(io.StringIO(csv_data))

        # Filter for chemical synapses ('S') originating from the 'ASJ' neuron
        asj_outputs = df[(df['neuron_1'] == 'ASJ') & (df['type'] == 'S')]

        if asj_outputs.empty:
            print("No chemical synapse outputs found for neuron ASJ.")
            return

        # Group by the postsynaptic neuron ('neuron_2') and sum the number of synapses ('N')
        synapse_counts = asj_outputs.groupby('neuron_2')['N'].sum().sort_values(ascending=False)

        # Get the main target (the one with the most synapses)
        main_target_cell = synapse_counts.idxmax()
        max_synapses = synapse_counts.max()

        print("\n--- Analysis Results ---")
        print(f"The main projection target of ASJ axons is: {main_target_cell}")
        print(f"Total synapses from ASJ to {main_target_cell}: {max_synapses}")
        
        print("\nBreakdown of synapse counts for all ASJ targets:")
        # To satisfy "output each number in the final equation"
        for target, count in synapse_counts.items():
            print(f"ASJ -> {target}: {count} synapses")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred while fetching the data: {e}")
    except Exception as e:
        print(f"An error occurred during data processing: {e}")

if __name__ == '__main__':
    find_asj_main_target()