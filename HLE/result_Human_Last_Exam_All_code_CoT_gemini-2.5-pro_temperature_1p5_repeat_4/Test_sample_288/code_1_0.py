import pandas as pd
# It may be necessary to install the library first:
# pip install c-elegans-wiring
from c_elegans_wiring.wiring_data import get_dataset

def find_main_asj_target():
    """
    Finds the main synaptic target of ASJ neurons in the adult C. elegans.
    """
    try:
        # 1. Load the adult hermaphrodite connectome dataset
        connections_df = get_dataset()

        # 2. Isolate connections originating from ASJ neurons (ASJL and ASJR)
        asj_outputs = connections_df[connections_df['pre'].isin(['ASJL', 'ASJR'])]

        # 3. Group by postsynaptic neuron and sum the number of synapses
        synapse_counts_by_target = asj_outputs.groupby('post')['N'].sum()

        if synapse_counts_by_target.empty:
            print("No synaptic targets found for ASJ neurons in the dataset.")
            return

        # 4. Identify the main target neuron with the maximum number of synapses
        main_target_cell = synapse_counts_by_target.idxmax()
        total_synapses = synapse_counts_by_target.max()

        # 5. Get the individual synapse counts from ASJL and ASJR to the main target
        asjl_to_target = asj_outputs[(asj_outputs['pre'] == 'ASJL') & (asj_outputs['post'] == main_target_cell)]['N'].sum()
        asjr_to_target = asj_outputs[(asj_outputs['pre'] == 'ASJR') & (asj_outputs['post'] == main_target_cell)]['N'].sum()

        # 6. Print the results
        print(f"Finding the main projection target of ASJ axons by synapse number...")
        print(f"Synapses from ASJL to {main_target_cell}: {asjl_to_target}")
        print(f"Synapses from ASJR to {main_target_cell}: {asjr_to_target}")
        print("-" * 30)
        print(f"Total synapses from ASJ to {main_target_cell}: {asjl_to_target} + {asjr_to_target} = {total_synapses}")
        print("\nConclusion:")
        print(f"The main projection target of ASJ axons is the cell '{main_target_cell}' with a total of {total_synapses} synapses.")

    except ImportError:
        print("Please install the required library by running: pip install c-elegans-wiring")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    find_main_asj_target()
    print("\n<<<AIY>>>")