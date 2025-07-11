import pandas as pd
import re

def solve_asj_target():
    """
    This function identifies the main projection target of ASJ axons
    in the adult C. elegans connectome by synapse number.
    """
    try:
        from c_elegans_wiring.synapses import query_synapses
    except ImportError:
        print("Please install the required library first by running: pip install c.elegans-wiring pandas")
        return

    # Load the synapse data from the Witvliet et al. 2021 dataset.
    synapses_df = query_synapses()

    # Define the presynaptic neurons of interest: ASJ Left and ASJ Right.
    asj_presynaptic_neurons = ['ASJL', 'ASJR']

    # Filter the dataframe to get all connections originating from ASJ neurons.
    asj_outputs = synapses_df[synapses_df['pre_cell'].isin(asj_presynaptic_neurons)].copy()

    # Define a function to get the cell "type" from a full cell name by removing L/R suffixes.
    def get_cell_type(cell_name):
        return re.sub(r'[LR]$', '', str(cell_name))

    # Apply this function to create a new column for the target cell type.
    asj_outputs['post_cell_type'] = asj_outputs['post_cell'].apply(get_cell_type)

    # Group by the target cell type and sum the number of synapses ('N').
    target_counts = asj_outputs.groupby('post_cell_type')['N'].sum().reset_index()

    # Sort the results to find the target with the highest number of synapses.
    sorted_targets = target_counts.sort_values(by='N', ascending=False)

    if sorted_targets.empty:
        print("Could not find any projection targets for ASJ neurons.")
        return

    # Get the top target's information.
    main_target_row = sorted_targets.iloc[0]
    main_target_type = main_target_row['post_cell_type']

    # Find the individual synaptic connections from ASJ to this main target type.
    main_target_details = asj_outputs[asj_outputs['post_cell_type'] == main_target_type]
    
    print(f"To find the main target of ASJ, we sum the synapses to each target cell type.")
    print(f"The top target cell type is {main_target_type}.")
    print("\nThe equation for the total number of synapses from ASJ to this target is:")

    equation_parts = []
    total_synapses = 0
    # Sort details for a consistent output order
    main_target_details_sorted = main_target_details.sort_values(by=['pre_cell', 'post_cell'])

    for _, row in main_target_details_sorted.iterrows():
        pre_cell = row['pre_cell']
        post_cell = row['post_cell']
        count = row['N']
        # Each number in the "equation" is a specific connection's synapse count
        equation_parts.append(str(count))
        total_synapses += count
        print(f"Synapses from {pre_cell} to {post_cell}: {count}")

    equation_str = " + ".join(equation_parts)
    print(f"\nFinal calculation: {equation_str} = {total_synapses}")
    print("\n---")
    print(f"Conclusion: The main projection target of ASJ axons is the {main_target_type} cell type, with a total of {total_synapses} synapses.")

if __name__ == '__main__':
    solve_asj_target()