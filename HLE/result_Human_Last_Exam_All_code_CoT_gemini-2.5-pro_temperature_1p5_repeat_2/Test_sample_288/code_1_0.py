import pandas as pd
try:
    import c_elegans_connectome as c
except ImportError:
    print("Please install the required library first by running: pip install c-elegans-connectome")
    exit()

def find_main_target():
    """
    Finds the main synaptic target of ASJ neurons in the adult C. elegans connectome.
    """
    # Load the connectome data for the adult hermaphrodite
    print("Loading C. elegans connectome data...")
    con = c.Connectome()

    # Define the presynaptic neurons of interest (ASJ left and right)
    source_neurons = ['ASJL', 'ASJR']

    # Get the connectivity matrix, which includes all connections
    # This returns a pandas DataFrame
    all_connections = con.get_connectivity()

    # Filter for connections originating from ASJ neurons
    asj_outputs = all_connections[all_connections['pre_cell'].isin(source_neurons)]

    # Group by the target neuron ('post_cell') and sum the number of synapses ('N')
    # This aggregates the synapses from both ASJL and ASJR to each target cell
    target_synapse_counts = asj_outputs.groupby('post_cell')['N'].sum()

    # Find the target cell with the maximum number of synapses
    main_target_cell = target_synapse_counts.idxmax()
    max_synapses = target_synapse_counts.max()

    # Display the result
    print(f"\nAnalysis of ASJ ('{source_neurons[0]}' and '{source_neurons[1]}') axonal targets:")
    print(f"Total synapses from ASJL to {main_target_cell}: {asj_outputs[(asj_outputs['pre_cell'] == 'ASJL') & (asj_outputs['post_cell'] == main_target_cell)]['N'].sum()}")
    print(f"Total synapses from ASJR to {main_target_cell}: {asj_outputs[(asj_outputs['pre_cell'] == 'ASJR') & (asj_outputs['post_cell'] == main_target_cell)]['N'].sum()}")
    print(f"\nThe main synaptic target of ASJ is the neuron pair {main_target_cell}, receiving a total of {int(max_synapses)} synapses.")
    print(f"Final Equation: ASJ -> {main_target_cell} ({int(max_synapses)} synapses)")


if __name__ == "__main__":
    find_main_target()

<<<AIA>>>