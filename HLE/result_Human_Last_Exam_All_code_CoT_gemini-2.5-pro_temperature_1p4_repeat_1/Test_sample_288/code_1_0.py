import pandas as pd
import io
import requests

def get_neuron_class(name):
    """Removes a terminal 'L' or 'R' from a neuron name to get its class."""
    if isinstance(name, str):
        if name.endswith('L') or name.endswith('R'):
            # Check if the character before L/R is a letter to avoid changing names like 'CANL' to 'CAN' if 'CAN' is a different cell.
            # In C. elegans nomenclature, this is generally safe.
            return name[:-1]
    return name

def find_main_asj_target():
    """
    Finds the main projection target of ASJ axons by synapse number
    in the C. elegans connectome.
    """
    # URL for the connectome data from Witvliet et al., 2021
    url = "https://wormwiring.org/db/latest/all_connections.csv"
    print("Fetching C. elegans connectome data...")

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()

        # Load the data into a pandas DataFrame
        col_names = ['pre', 'post', 'type', 'N', 'neurotransmitter']
        df = pd.read_csv(io.StringIO(response.text), header=0, names=col_names)

        # Filter for chemical synapses (type 'S') from ASJ neurons ('ASJL', 'ASJR')
        asj_neurons = ['ASJL', 'ASJR']
        asj_outputs = df[(df['pre'].isin(asj_neurons)) & (df['type'] == 'S')].copy()

        if asj_outputs.empty:
            print("No synaptic outputs found for ASJ neurons.")
            return

        # Create a new column for the postsynaptic neuron class
        asj_outputs['post_class'] = asj_outputs['post'].apply(get_neuron_class)

        # Group by the postsynaptic neuron class and sum the number of synapses
        class_synapse_counts = asj_outputs.groupby('post_class')['N'].sum()

        # Find the class with the maximum number of synapses
        main_target_class = class_synapse_counts.idxmax()
        total_synapses = class_synapse_counts.max()

        # Get the detailed connections to the main target class
        main_target_details = asj_outputs[asj_outputs['post_class'] == main_target_class]
        
        # Get synapse counts from ASJL and ASJR to their respective partners in the target class
        synapses_from_asjl = main_target_details[main_target_details['pre'] == 'ASJL']['N'].sum()
        synapses_from_asjr = main_target_details[main_target_details['pre'] == 'ASJR']['N'].sum()
        
        # Print the final result and the equation
        print("\n--- Analysis Complete ---")
        print(f"The main projection target of ASJ axons is the {main_target_class} neuron class.")
        print("The total number of synapses is calculated as follows:")
        print(f"Synapses from ASJL to {main_target_class}L: {synapses_from_asjl}")
        print(f"Synapses from ASJR to {main_target_class}R: {synapses_from_asjr}")
        print("\nFinal Equation:")
        print(f"{synapses_from_asjl} + {synapses_from_asjr} = {total_synapses}")


    except requests.exceptions.RequestException as e:
        print(f"Failed to retrieve data from {url}. Error: {e}")
    except Exception as e:
        print(f"An error occurred during the analysis: {e}")

if __name__ == '__main__':
    find_main_asj_target()