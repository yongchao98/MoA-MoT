# First, ensure you have the necessary libraries installed:
# pip install cbi_python pandas

from cbi_python import connectome
import pandas as pd

try:
    # 1. Load the C. elegans connectome data from Witvliet et al., 2021
    c = connectome.Connectome()
    # By default, this loads the 'witvliet_2021' dataset.
    # We get the data as a pandas DataFrame.
    df = c.get_connectivity_df()

    # 2. Filter for connections where ASJ is the presynaptic neuron
    # The 'pre_cell' column contains the sending neuron.
    asj_outputs = df[df['pre_cell'] == 'ASJ']

    # 3. Group by the postsynaptic target and sum the synapse counts
    # The 'post_cell' column contains the receiving neuron, and 'N' is the synapse count.
    synapse_counts_by_target = asj_outputs.groupby('post_cell')['N'].sum()

    # 4. Find the target with the maximum number of synapses
    if not synapse_counts_by_target.empty:
        main_target_cell = synapse_counts_by_target.idxmax()
        max_synapses = synapse_counts_by_target.max()

        # 5. Print the final result, showing each component of the finding.
        print(f"Analysis of connections from the presynaptic neuron: ASJ")
        print(f"The main postsynaptic target cell is: {main_target_cell}")
        print(f"The total number of synapses to this cell is: {int(max_synapses)}")
        print(f"\nTherefore, the main projection target of ASJ axons is {main_target_cell} with {int(max_synapses)} synapses.")
    else:
        print("Could not find any synaptic outputs from ASJ in the dataset.")

except ImportError:
    print("Please install the required libraries by running: pip install cbi_python pandas")
except Exception as e:
    print(f"An error occurred: {e}")
    print("This may be due to a network issue or a problem with the cbi_python library. Please check your internet connection and library installation.")
