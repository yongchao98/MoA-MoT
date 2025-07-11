import pandas as pd
from neuprint import Client

def find_asj_main_target():
    """
    Connects to the C. elegans neuPrint server to find the main synaptic
    target of the ASJ neuron pair.
    """
    try:
        # Step 1: Connect to the neuPrint server and select the C. elegans dataset
        # This public dataset does not require an authentication token.
        c = Client('https://neuprint.org', dataset='c_elegans_isbi_2021:v1')
        c.fetch_version()
        print("Successfully connected to neuPrint server.")
    except Exception as e:
        print(f"Error connecting to neuPrint: {e}")
        print("Please check your internet connection and the server status.")
        return

    # Step 2: Define the presynaptic neurons and the Cypher query
    pre_neurons = ['ASJL', 'ASJR']
    cypher_query = f"""
    MATCH (pre:Neuron)-[r:ConnectsTo]->(post:Neuron)
    WHERE pre.name IN {pre_neurons}
    RETURN pre.name AS presynaptic, post.name AS postsynaptic, r.weight AS synapses
    """

    try:
        # Step 3: Fetch the connection data
        results_df = c.fetch_cypher(cypher_query)
        if results_df.empty:
            print(f"No connections found for neurons: {pre_neurons}")
            return

        # Step 4: Aggregate the results to find the total synapses for each target
        total_synapses = results_df.groupby('postsynaptic')['synapses'].sum().sort_values(ascending=False).astype(int)

        # Step 5: Identify the main target
        main_target = total_synapses.index[0]
        total_count = total_synapses.iloc[0]

        # Get the individual synapse counts from ASJL and ASJR for the main target
        main_target_connections = results_df[results_df['postsynaptic'] == main_target]
        asjl_row = main_target_connections[main_target_connections['presynaptic'] == 'ASJL']
        asjr_row = main_target_connections[main_target_connections['presynaptic'] == 'ASJR']
        
        asjl_count = asjl_row['synapses'].iloc[0] if not asjl_row.empty else 0
        asjr_count = asjr_row['synapses'].iloc[0] if not asjr_row.empty else 0

        # Step 6: Print the results clearly
        print("\n---")
        print(f"The main projection target of ASJ in terms of synapse number is: {main_target}")
        print("\nCalculation:")
        # Print the final equation with each number
        print(f"ASJL -> {main_target} ({asjl_count}) + ASJR -> {main_target} ({asjr_count}) = {total_count} synapses")
        print("\nFull list of targets ranked by synapse number:")
        print(total_synapses.to_string())

    except Exception as e:
        print(f"An error occurred during data processing: {e}")

if __name__ == "__main__":
    find_asj_main_target()
<<<AVG>>>