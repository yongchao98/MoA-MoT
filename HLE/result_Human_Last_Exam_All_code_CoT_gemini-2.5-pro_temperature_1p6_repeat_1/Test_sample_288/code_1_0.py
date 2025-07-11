import requests
import pandas as pd
import re

def find_asj_main_target():
    """
    Queries the wormwiring.org connectome to find the main synaptic target 
    cell type of ASJ neurons based on synapse number.
    """
    # Step 1: Define the API endpoint and the Cypher query.
    api_url = "https://wormwiring.org/db/api/query"
    # This query retrieves all presynaptic connections from ASJ neurons.
    cypher_query = {
        "cypher": "MATCH (n1:Neuron)-[r:Synapse]->(n2:Neuron) WHERE n1.name IN ['ASJL', 'ASJR'] RETURN n1.name AS Presynaptic, n2.name AS Postsynaptic, r.N AS Synapses"
    }

    print("Querying the C. elegans connectome database (wormwiring.org)...")

    try:
        # Step 2: Make the POST request to the API.
        response = requests.post(api_url, json=cypher_query)
        response.raise_for_status()  # Raise an exception for HTTP errors
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not retrieve data from the API. {e}")
        return

    try:
        data = response.json()
    except requests.exceptions.JSONDecodeError:
        print("Error: Could not decode JSON from the API response.")
        return

    if 'data' not in data or not data['data']:
        print("Query returned no data. No synaptic targets found for ASJ.")
        return

    # Create a pandas DataFrame for easy data manipulation.
    df = pd.DataFrame(data['data'], columns=data['columns'])
    # Ensure the 'Synapses' column is numeric for calculations.
    df['Synapses'] = pd.to_numeric(df['Synapses'])

    # Step 3: Group neurons by cell type (e.g., RIAL and RIAR become RIA).
    # This is done by removing the 'L' or 'R' suffix from the neuron name.
    df['Postsynaptic_Type'] = df['Postsynaptic'].apply(lambda name: re.sub(r'(L|R)$', '', name))

    # Step 4: Calculate the total number of synapses for each postsynaptic cell type.
    type_sums = df.groupby('Postsynaptic_Type')['Synapses'].sum().sort_values(ascending=False)

    if type_sums.empty:
        print("No synaptic targets found for ASJ.")
        return

    # Step 5: Identify the main target cell type and the total synapse count.
    main_target_type = type_sums.index[0]
    total_synapses = type_sums.iloc[0]

    # Get the detailed breakdown for the main target type.
    main_target_df = df[df['Postsynaptic_Type'] == main_target_type].sort_values(by='Postsynaptic')

    # Step 6: Display the final results.
    print(f"\nThe main projection target of ASJ axons is the {main_target_type} cell type.")
    print("\nBreakdown of the synaptic connections:")

    equation_parts = []
    for _, row in main_target_df.iterrows():
        synapse_count = row['Synapses']
        print(f"- From {row['Presynaptic']} to {row['Postsynaptic']}: {synapse_count} synapses")
        equation_parts.append(str(synapse_count))

    # Print the final calculation as a clear equation.
    equation_str = f"{' + '.join(equation_parts)} = {total_synapses}"
    print(f"\nThe total number of synapses from ASJ to {main_target_type} is calculated as:")
    print(equation_str)

if __name__ == '__main__':
    find_asj_main_target()