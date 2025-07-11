import celegans_connectome as cc
import pandas as pd

try:
    # 1. Load the connectome data from Witvliet et al., 2021
    c = cc.Connectome()
    chem_adj = c.get_adjacency_matrix(kind='chemical')

    # 2. Define the source neurons
    source_neurons = ['ASJL', 'ASJR']

    # 3. Sum the outputs from the ASJ neuron pair to all other neurons
    # We use .reindex to handle cases where a neuron might be missing from the data
    # and fill with 0.
    asj_outputs = chem_adj.reindex(source_neurons).fillna(0).sum(axis=0)

    # 4. Find the neuron that receives the most synapses from ASJ
    main_target_cell = asj_outputs.idxmax()
    total_synapses = int(asj_outputs.max())

    # Get the individual synapse counts for the final equation
    synapses_from_asjl = int(chem_adj.loc['ASJL', main_target_cell])
    synapses_from_asjr = int(chem_adj.loc['ASJR', main_target_cell])

    # 5. Print the results and the equation
    print(f"The main projection target of ASJ axons is the neuron: {main_target_cell}")
    print("\nBreakdown of synapse counts:")
    print(f"ASJL -> {main_target_cell}: {synapses_from_asjl}")
    print(f"ASJR -> {main_target_cell}: {synapses_from_asjr}")
    print("\nFinal Equation:")
    print(f"Total ASJ -> {main_target_cell} = {synapses_from_asjl} + {synapses_from_asjr} = {total_synapses}")

except Exception as e:
    print(f"An error occurred. Please ensure you have the 'celegans_connectome' library installed (`pip install celegans_connectome`).")
    print(f"Error details: {e}")
