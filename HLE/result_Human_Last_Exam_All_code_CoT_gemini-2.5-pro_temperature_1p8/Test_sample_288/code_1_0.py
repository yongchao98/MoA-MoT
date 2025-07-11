import pandas as pd
import io
import requests

def find_main_asj_target():
    """
    This script finds the main synaptic target of the ASJ neurons in C. elegans.
    It fetches the connectome data, filters for chemical synapses from ASJL and ASJR,
    aggregates the synapse counts for each target, and identifies the target with the highest count.
    """
    # URL for the C. elegans hermaphrodite full connectome edgelist
    url = "http://www.openconnecto.me/data/celegans/connectome/herm_full_edgelist.csv"
    
    try:
        # Step 1: Fetch the data
        print("Fetching C. elegans connectome data...")
        response = requests.get(url)
        response.raise_for_status()  # This will raise an HTTPError for bad responses (4xx or 5xx)
        data = response.content.decode('utf-8')
        
        # Load the data into a pandas DataFrame. The columns are 'source', 'target', 'type', 'weight'
        df = pd.read_csv(io.StringIO(data))

        # Step 2: Filter for chemical synapses ('S') from ASJ neurons
        asj_neurons = ['ASJL', 'ASJR']
        asj_projections = df[
            (df['source'].isin(asj_neurons)) & 
            (df['type'] == 'S')
        ]

        if asj_projections.empty:
            print("No chemical synapse projections found for ASJ neurons in the dataset.")
            return

        # Step 3: Group by target cell and sum the synapse weights
        # This calculates the total number of synapses from the ASJ pair to each target.
        synapse_counts_per_target = asj_projections.groupby('target')['weight'].sum().sort_values(ascending=False)

        # Step 4: Identify the main target (cell with the most synapses)
        if synapse_counts_per_target.empty:
            print("Could not aggregate synapse counts.")
            return

        main_target_cell = synapse_counts_per_target.index[0]
        max_synapses = synapse_counts_per_target.iloc[0]

        # Step 5: Print the final result, including the numbers used
        print(f"\nAnalysis complete.")
        print(f"The main projection target of ASJ axons is the cell: {main_target_cell}")
        print(f"Total number of synapses from ASJ (ASJL & ASJR) to {main_target_cell}: {max_synapses}")
        
        # Note on tie-breaking:
        # Check if there is a tie for the top spot
        top_targets = synapse_counts_per_target[synapse_counts_per_target == max_synapses]
        if len(top_targets) > 1:
            print("\nNote: There is a tie. Other cells with the same number of synapses include:")
            for cell in top_targets.index[1:]:
                print(f"- {cell} ({top_targets[cell]} synapses)")


    except requests.exceptions.RequestException as e:
        print(f"Error: Failed to retrieve data from {url}.")
        print(f"Please check your internet connection or the URL's validity. Details: {e}")
    except Exception as e:
        print(f"An unexpected error occurred during data processing: {e}")

if __name__ == "__main__":
    find_main_asj_target()