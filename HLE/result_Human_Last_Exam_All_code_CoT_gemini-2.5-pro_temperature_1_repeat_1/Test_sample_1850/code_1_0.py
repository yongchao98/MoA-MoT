def count_saints_in_paradise_lost():
    """
    This script identifies and counts the historical saints mentioned by name
    or clear allusion in John Milton's "Paradise Lost".
    """
    
    # Based on analysis of the text, particularly Book III, lines 473-477,
    # two historical saints are mentioned in a critique of monastic orders.
    saints = {
        "St. Dominic": "Mentioned by name: '...weeds of DOMINIC...'",
        "St. Francis": "Referenced by his order: '...in Franciscan think to pass disguis'd...'"
    }
    
    saint_names = list(saints.keys())
    count = len(saint_names)
    
    print("The historical saints mentioned in Milton's Paradise Lost are:")
    for name, reference in saints.items():
        print(f"- {name} ({reference})")
        
    print("\nFinal Calculation:")
    
    # To fulfill the request of showing the equation, we list the names.
    # We extract the core name (e.g., "Dominic" from "St. Dominic").
    equation_parts = [name.split(' ')[1] for name in saint_names]
    equation_str = " + ".join(equation_parts)
    
    print(f"{equation_str} = {count}")

count_saints_in_paradise_lost()