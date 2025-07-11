import pandas as pd
import io
import requests

def find_asj_main_target():
    """
    This function downloads the C. elegans connectome data, finds the main
    synaptic target of the ASJ neurons, and prints the result.
    """
    # Step 1: Fetch the connectome data from wormwiring.org
    url = "http://wormwiring.org/db/dataset/1/download/N2U_neuron_connections.csv"
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        csv_data = response.content
    except requests.exceptions.RequestException as e:
        print(f"Failed to download data: {e}")
        return

    # Use pandas to read the CSV data
    df = pd.read_csv(io.StringIO(csv_data.decode('utf-8')))
    
    # We are interested in chemical synapses, which are of type 'C'
    df_chem = df[df['N2U_type'] == 'C'].copy()

    # Step 2: Filter for connections originating from ASJ neurons (ASJL and ASJR)
    source_neurons = ['ASJL', 'ASJR']
    asj_projections = df_chem[df_chem['pre'].isin(source_neurons)]

    if asj_projections.empty:
        print("No connections found for ASJ neurons in the dataset.")
        return
        
    # Step 3: Aggregate synapse counts for each target neuron
    # We group by the postsynaptic neuron ('post') and sum the synapse weights ('N2U_weight')
    target_synapse_counts = asj_projections.groupby('post')['N2U_weight'].sum()

    # Step 4: Identify the main target neuron with the maximum number of synapses
    main_target_neuron = target_synapse_counts.idxmax()
    total_synapses_to_target = target_synapse_counts.max()

    # Step 5: Get the breakdown from ASJL and ASJR to the main target
    asjl_to_target_synapses = asj_projections[
        (asj_projections['pre'] == 'ASJL') & (asj_projections['post'] == main_target_neuron)
    ]['N2U_weight'].sum()

    asjr_to_target_synapses = asj_projections[
        (asj_projections['pre'] == 'ASJR') & (asj_projections['post'] == main_target_neuron)
    ]['N2U_weight'].sum()

    # Print the final result as an equation
    print(f"The main projection target of ASJ is {main_target_neuron}.")
    print("\nThis is based on the total number of chemical synapses received from ASJL and ASJR:")
    print(f"\nTarget: {main_target_neuron} ({total_synapses_to_target} total synapses)")
    print("Equation:")
    print(f"{asjl_to_target_synapses} (from ASJL) + {asjr_to_target_synapses} (from ASJR) = {total_synapses_to_target}")

if __name__ == "__main__":
    find_asj_main_target()