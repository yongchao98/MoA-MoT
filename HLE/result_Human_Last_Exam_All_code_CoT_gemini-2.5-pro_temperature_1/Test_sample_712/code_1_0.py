def find_triboluminescent_bacterium():
    """
    This function simulates a query to a knowledge base to find the genus
    of bacterium known for producing a triboluminescent compound.
    """
    # A simplified knowledge base representation.
    # Data is structured to link phenomena to organisms and their products.
    knowledge_base = {
        'phenomena': {
            'triboluminescence': {
                'source_organism': 'bacterium',
                'organism_genus': 'Streptomyces',
                'compound': 'avermectin'
            }
        }
    }

    # Retrieve the genus from the knowledge base.
    genus = knowledge_base['phenomena']['triboluminescence']['organism_genus']
    compound = knowledge_base['phenomena']['triboluminescence']['compound']

    print(f"The genus of bacterium that produces a triboluminescent compound is: {genus}")
    print(f"The specific compound is {compound}, which exhibits triboluminescence when its crystals are crushed.")

# Execute the function to find and print the answer.
find_triboluminescent_bacterium()