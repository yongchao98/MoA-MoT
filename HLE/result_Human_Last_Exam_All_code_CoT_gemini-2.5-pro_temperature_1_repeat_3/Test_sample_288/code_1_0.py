import pandas as pd
try:
    # Use the cenelegans library to fetch the data
    from cenelegans.data import get_data
    # Load the chemical synapse data for the adult hermaphrodite
    chem_df = get_data('herm_chem')
except (ImportError, ModuleNotFoundError):
    print("cenelegans library not found. Falling back to a cached remote URL.")
    print("For a more reliable and up-to-date connection, please run: pip install cenelegans")
    # Fallback to a URL if cenelegans is not installed
    url = 'https://raw.githubusercontent.com/openconnectome/data_wrangler/master/wrangler/data/Celegans/hermaphrodite/chem.csv'
    chem_df = pd.read_csv(url)
    # Rename columns to match the cenelegans library for consistency
    chem_df = chem_df.rename(columns={'Source': 'pre', 'Target': 'post', 'Weight': 'N'})

# Define the source neurons (ASJ left and right)
source_neurons = ['ASJL', 'ASJR']

# Filter for connections originating from ASJ neurons
asj_connections = chem_df[chem_df['pre'].isin(source_neurons)]

# Group by the target (postsynaptic) neuron and sum the number of synapses ('N')
target_totals = asj_connections.groupby('post')['N'].sum()

# Find the target neuron with the maximum number of synapses
if not target_totals.empty:
    main_target_cell = target_totals.idxmax()
    total_synapses = target_totals.max()

    # Get the synapse count from ASJL to the main target
    asjl_to_target_synapses = asj_connections[
        (asj_connections['pre'] == 'ASJL') & (asj_connections['post'] == main_target_cell)
    ]['N'].sum()

    # Get the synapse count from ASJR to the main target
    asjr_to_target_synapses = asj_connections[
        (asj_connections['pre'] == 'ASJR') & (asj_connections['post'] == main_target_cell)
    ]['N'].sum()

    print(f"The main projection target of ASJ axons is {main_target_cell} with a total of {total_synapses} synapses.")
    print("\nBreakdown:")
    print(f"ASJL -> {main_target_cell}: {asjl_to_target_synapses} synapses")
    print(f"ASJR -> {main_target_cell}: {asjr_to_target_synapses} synapses")
    print(f"Total = {asjl_to_target_synapses} + {asjr_to_target_synapses} = {total_synapses} synapses")
else:
    print("Could not find any connections from ASJ neurons in the dataset.")

<<<PVQ>>>