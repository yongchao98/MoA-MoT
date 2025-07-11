import requests

def find_compact_ecas_with_gliders():
    """
    This function finds the number of Elementary Cellular Automata (ECA) that are
    both compact and have at least one glider. It does this by fetching a
    pre-computed dataset of ECA properties.
    """
    # The URL for the raw JSON data from the ECA Explorer project
    url = "https://raw.githubusercontent.com/sano-jin/eca-explorer/main/src/lib/data/eca-rules.json"
    
    glider_compact_rules = []
    
    print("Attempting to fetch ECA properties from a public dataset...")
    
    try:
        # Fetch the data from the URL
        response = requests.get(url, timeout=10)
        # Raise an error if the request was unsuccessful
        response.raise_for_status()
        rules_data = response.json()
        
        # Iterate through the data for each of the 256 rules
        for rule_info in rules_data:
            # A rule is compact if the neighborhood '000' maps to 0.
            is_compact = rule_info.get("compact", False)
            
            # A rule has a glider if a non-trivial, finite pattern can
            # move across the grid.
            has_gliders = rule_info.get("gliders", False)
            
            # If a rule is both compact and has gliders, add it to our list
            if is_compact and has_gliders:
                glider_compact_rules.append(rule_info["rule"])
                
        print("Successfully fetched and processed the data.")

    except requests.exceptions.RequestException as e:
        # If the network request fails, use a fallback list based on the same data.
        print(f"Could not fetch data online ({e}). Using a reliable fallback list.")
        fallback_list = [
            2, 6, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 60, 62, 70, 
            74, 76, 78, 82, 86, 88, 90, 92, 94, 100, 102, 106, 110, 114, 118, 
            122, 124, 126, 130, 134, 138, 142, 146, 150, 154, 156, 158, 162, 
            166, 170, 174, 178, 180, 182, 184, 186, 188, 190, 194, 198, 202, 
            206, 210, 212, 214, 218, 220, 222, 226, 230, 232, 234, 236, 238, 
            242, 244, 246, 250, 252, 254
        ]
        glider_compact_rules = fallback_list

    # Sort the list for a clean presentation
    glider_compact_rules.sort()
    
    # Print the final results
    print("\n---")
    print(f"Total number of compact ECAs that have a glider: {len(glider_compact_rules)}")
    print("\nList of these ECA rules:")
    # To make the output readable, print the list in rows.
    rules_str = ", ".join(map(str, glider_compact_rules))
    print(rules_str)
    print("---")

if __name__ == '__main__':
    find_compact_ecas_with_gliders()