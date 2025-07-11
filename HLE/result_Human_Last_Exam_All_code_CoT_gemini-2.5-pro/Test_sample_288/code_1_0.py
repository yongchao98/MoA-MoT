# First, ensure you have the necessary libraries installed:
# pip install celegans-connectome pandas

import pandas as pd
import celegans_connectome as ce

def find_asj_main_target():
    """
    This script analyzes the C. elegans connectome to find the main synaptic
    target of the ASJ neuron class.
    """
    try:
        # Load the neuron-to-neuron synapse data from the Witvliet et al. 2021 dataset.
        # This dataset is a comprehensive connectome of an adult C. elegans hermaphrodite.
        print("Loading C. elegans connectome data (Witvliet et al., 2021)...")
        connectome_df = ce.get_connectome_data('witvliet_2021', 'neuron_to_neuron_synapses')
        print("Data loaded successfully.")

        # Filter the data to get all connections originating from the ASJ neuron pair (ASJL and ASJR).
        asj_outputs = connectome_df[connectome_df['pre_cell'].isin(['ASJL', 'ASJR'])].copy()

        # Define a function to determine the neuron class from a specific neuron name
        # (e.g., 'AIYL' -> 'AIY', 'AVAL' -> 'AVA').
        def get_cell_class(cell_name):
            if cell_name.endswith('L') or cell_name.endswith('R'):
                # For paired neurons, the class is the name without the L/R suffix.
                return cell_name[:-1]
            # For unpaired neurons, the class is the name itself.
            return cell_name

        # Create a new column for the postsynaptic cell class to group the targets.
        asj_outputs['post_cell_class'] = asj_outputs['post_cell'].apply(get_cell_class)

        # Group by the target cell class and sum the synapse counts ('N').
        target_class_synapses = asj_outputs.groupby('post_cell_class')['N'].sum()

        # Find the target class that receives the maximum number of synapses.
        main_target_class = target_class_synapses.idxmax()
        max_synapses = target_class_synapses.max()

        print("\n--------------------------------------------------")
        print(f"The main projection target of ASJ axons is the '{main_target_class}' neuron class.")
        print("--------------------------------------------------\n")

        # To show the final equation, we find the individual connections to this target class.
        print(f"This is based on the following synaptic connections from ASJ to the {main_target_class} class:")
        
        # Filter for connections to the main target class
        main_target_connections = asj_outputs[asj_outputs['post_cell_class'] == main_target_class]
        
        # Get the individual connection counts
        connection_counts = main_target_connections.groupby(['pre_cell', 'post_cell'])['N'].sum()
        
        equation_parts = []
        for (pre_cell, post_cell), count in connection_counts.items():
            print(f"- Synapses from {pre_cell} to {post_cell}: {count}")
            equation_parts.append(str(count))

        # Format and print the final equation.
        equation_str = " + ".join(equation_parts)
        print(f"\nThe final equation for the total synapse count is:")
        print(f"{equation_str} = {max_synapses}")


    except Exception as e:
        print("\nAn error occurred. Please ensure you have the required libraries installed.")
        print("You can install them by running this command in your terminal:")
        print("pip install celegans-connectome pandas")
        print(f"\nError details: {e}")

if __name__ == "__main__":
    find_asj_main_target()
<<<AIY>>>