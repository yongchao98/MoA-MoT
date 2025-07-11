import pandas as pd
import io
import requests

def find_main_synaptic_target():
    """
    Finds the main synaptic target of a given source neuron pair in the C. elegans connectome.
    """
    # The URL points to a CSV file of the adult C. elegans chemical synapse connectome.
    url = "https://raw.githubusercontent.com/connectome-synapse/c-elegans-connectome-data/main/data/N2U_neuron_connection.csv"
    
    print(f"Fetching connectome data from {url}...")
    
    try:
        # Download the data
        response = requests.get(url)
        # Raise an error if the download fails
        response.raise_for_status()
        csv_data = response.text
        print("Data fetched successfully.")
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not fetch data from the URL. {e}")
        return

    # Load the data into a pandas DataFrame
    # The columns are 'origin' (presynaptic), 'target' (postsynaptic), and 'number' (synapse count)
    df = pd.read_csv(io.StringIO(csv_data))
    
    # The source neurons are the left/right pair of ASJ
    source_neurons = ['ASJL', 'ASJR']
    print(f"\nAnalyzing projections from source neurons: {', '.join(source_neurons)}")

    # Filter for all connections originating from ASJ neurons
    asj_projections = df[df['origin'].isin(source_neurons)]

    if asj_projections.empty:
        print("No projections found for ASJ neurons in the dataset.")
        return

    # Group connections by the target neuron and sum the synapse counts
    synapse_counts_per_target = asj_projections.groupby('target')['number'].sum()

    # Find the target neuron that receives the most synapses from ASJ
    main_target_neuron = synapse_counts_per_target.idxmax()
    max_synapses = synapse_counts_per_target.max()

    # Get the individual synapse counts from ASJL and ASJR to the main target
    synapses_from_asjl_series = asj_projections[(asj_projections['origin'] == 'ASJL') & (asj_projections['target'] == main_target_neuron)]['number']
    synapses_from_asjl = synapses_from_asjl_series.sum() if not synapses_from_asjl_series.empty else 0
    
    synapses_from_asjr_series = asj_projections[(asj_projections['origin'] == 'ASJR') & (asj_projections['target'] == main_target_neuron)]['number']
    synapses_from_asjr = synapses_from_asjr_series.sum() if not synapses_from_asjr_series.empty else 0

    print("\n--- Result ---")
    print(f"The main projection target of ASJ axons is the neuron '{main_target_neuron}'.")
    print(f"It receives a total of {max_synapses} synapses from the ASJ pair.")
    print(f"Equation: {synapses_from_asjl} (from ASJL) + {synapses_from_asjr} (from ASJR) = {max_synapses}")

if __name__ == "__main__":
    find_main_synaptic_target()