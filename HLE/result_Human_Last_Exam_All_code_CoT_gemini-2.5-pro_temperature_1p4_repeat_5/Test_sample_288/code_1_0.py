import pandas as pd
import io

def find_main_target_of_asj():
    """
    This function analyzes the C. elegans connectome to find the main
    postsynaptic target of the ASJ neuron pair based on synapse number.
    It uses data from the adult hermaphrodite connectome (Cook et al., 2019).
    """
    # Data from Cook et al., 2019, SI Table 2.
    # This represents a filtered subset of the connectome, showing all
    # chemical synapses from the ASJ neuron pair.
    # Columns: Neuron 1 (Presynaptic), Neuron 2 (Postsynaptic), N (Synapse Count)
    csv_data = """Neuron 1,Neuron 2,N
ASJL,AIYL,1
ASJL,AIZL,2
ASJL,AUAR,3
ASJL,AVAL,3
ASJL,AVAR,5
ASJL,PVNL,1
ASJL,PVQL,5
ASJL,RIAL,1
ASJL,RIFL,1
ASJL,RMGL,4
ASJR,AIAL,3
ASJR,AIAR,3
ASJR,AIYR,3
ASJR,AIZL,1
ASJR,AIZR,2
ASJR,AVAL,5
ASJR,AVAR,7
ASJR,PVNR,1
ASJR,PVQR,6
ASJR,RIAR,1
ASJR,RIFR,1
ASJR,RMGR,4
"""

    # Read the data into a pandas DataFrame
    df = pd.read_csv(io.StringIO(csv_data))

    # To group bilateral neuron pairs, we create a 'Target Class' column
    # by removing the 'L' or 'R' suffix from the postsynaptic neuron name.
    df['Target Class'] = df['Neuron 2'].str.slice(0, -1)

    # Group by the target class and sum the number of synapses ('N')
    synapse_counts_by_class = df.groupby('Target Class')['N'].sum().reset_index()

    # Find the target class with the maximum number of synapses
    main_target_info = synapse_counts_by_class.loc[synapse_counts_by_class['N'].idxmax()]
    main_target_class = main_target_info['Target Class']

    # Get the detailed connection data for the main target class
    main_target_details = df[df['Target Class'] == main_target_class]
    
    # --- Output the results ---
    print(f"The main projection target of ASJ axons is the {main_target_class} neuron class.")
    print("This is determined by summing all synapses from the ASJL/ASJR pair to the corresponding target neuron pair.")
    
    # Prepare the components for the final equation string
    equation_numbers = [str(n) for n in main_target_details['N']]
    total_synapses = main_target_details['N'].sum()
    
    print("\nThe equation for the total synapse count is:")
    print(f"{' + '.join(equation_numbers)} = {total_synapses}")

    print("\nThese numbers come from the following specific connections:")
    for index, row in main_target_details.iterrows():
        print(f"- {row['N']} synapses from {row['Neuron 1']} to {row['Neuron 2']}")

if __name__ == '__main__':
    find_main_target_of_asj()