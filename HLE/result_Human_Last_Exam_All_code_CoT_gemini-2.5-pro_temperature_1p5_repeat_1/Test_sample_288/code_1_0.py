import pandas as pd
import io
import requests

def find_main_target_neuron():
    """
    Fetches C. elegans connectome data to find the main synaptic target of ASJ neurons.
    """
    # The URL for the N2U adult hermaphrodite connectome data (chemical synapses)
    url = "http://www.wormwiring.org/db/DATA/N2U_neuron_connect.csv"
    print(f"Fetching connectome data from {url}...")

    try:
        response = requests.get(url)
        # Raise an HTTPError if the HTTP request returned an unsuccessful status code
        response.raise_for_status()
        
        # The data is embedded within <pre> tags in the HTML response. We need to extract it.
        # This makes the script robust to changes in the surrounding page structure.
        if '<pre>' in response.text:
            csv_content = response.text.split('<pre>')[1].split('</pre>')[0].strip()
        else:
            print("Could not find <pre> tag in the response. Assuming plain text.")
            csv_content = response.text

        # Use io.StringIO to treat the string data as a file for pandas
        # The first line of the actual data is a comment that pandas can skip.
        df = pd.read_csv(io.StringIO(csv_content), comment='#')

        # Clean up column names by removing leading/trailing whitespace
        df.columns = [col.strip() for col in df.columns]

        # Rename columns for easier access
        df = df.rename(columns={
            'Neuron 1': 'presynaptic_neuron',
            'Neuron 2': 'postsynaptic_neuron',
            'Nbr': 'number_of_synapses'
        })

        # Filter for connections where the presynaptic neuron is ASJL or ASJR
        asj_outputs = df[df['presynaptic_neuron'].isin(['ASJL', 'ASJR'])].copy()
        
        if asj_outputs.empty:
            print("No connections from ASJ neurons found in the dataset.")
            return

        print("\nFound connections from ASJ. Aggregating by target neuron...")

        # Group by the target neuron and sum the number of synapses
        target_synapse_counts = asj_outputs.groupby('postsynaptic_neuron')['number_of_synapses'].sum()

        # Find the neuron with the maximum number of synapses from ASJ
        main_target_cell = target_synapse_counts.idxmax()
        max_synapses = target_synapse_counts.max()

        print(f"\nAnalysis complete. The main target is {main_target_cell}.")
        print("-" * 40)
        print(f"The main projection target of ASJ axons is the cell '{main_target_cell}'.")
        
        # Breakdown of the synapse count from ASJL and ASJR to the main target
        asjl_synapses = asj_outputs[
            (asj_outputs['presynaptic_neuron'] == 'ASJL') & 
            (asj_outputs['postsynaptic_neuron'] == main_target_cell)
        ]['number_of_synapses'].sum()
        
        asjr_synapses = asj_outputs[
            (asj_outputs['presynaptic_neuron'] == 'ASJR') & 
            (asj_outputs['postsynaptic_neuron'] == main_target_cell)
        ]['number_of_synapses'].sum()

        print("\nThe total number of synapses from ASJ to this cell is calculated as follows:")
        print(f"Synapses from ASJL to {main_target_cell}: {asjl_synapses}")
        print(f"Synapses from ASJR to {main_target_cell}: {asjr_synapses}")
        print(f"Total synapses: {asjl_synapses} + {asjr_synapses} = {int(max_synapses)}")
        print("-" * 40)


    except requests.exceptions.RequestException as e:
        print(f"An error occurred while fetching the data: {e}")
    except (KeyError, IndexError, AttributeError) as e:
        print(f"An error occurred while processing the data. The data format may have changed. Error: {e}")

if __name__ == '__main__':
    find_main_target_neuron()