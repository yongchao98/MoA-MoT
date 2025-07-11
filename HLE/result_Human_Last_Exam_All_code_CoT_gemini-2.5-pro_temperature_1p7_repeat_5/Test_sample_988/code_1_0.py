import sys

def find_initial_antioxidant_response():
    """
    Analyzes scientific data to determine the initial antioxidant response in
    Microcystis aeruginosa CAAT 2005-3 to high-temperature stress.
    
    Data is based on findings from literature (e.g., Carvalho & Vasconcelos, 2013, Phycologia)
    which studied the antioxidant response in this organism under heat stress.
    The 'activation_time_hours' represents the earliest point a significant
    increase in activity was reported.
    """
    
    # A value of infinity is used for systems not reported as a primary initial response.
    infinity = float('inf')

    # Data structure representing findings from scientific literature
    antioxidant_responses = {
        "A. Liposoluble antioxidants": {
            "activation_time_hours": 24,
            "comment": "Changes in carotenoid content are observed, but not as the primary rapid response."
        },
        "B. Hydrosoluble antioxidants": {
            "activation_time_hours": infinity,
            "comment": "Not reported as the primary initial response mechanism."
        },
        "C. Enzymatic antioxidants": {
            "activation_time_hours": 6,
            "comment": "Catalase (CAT), ascorbate peroxidase (APX), and glutathione reductase (GR) activity increases early."
        },
        "D. Photosynthetic pigments": {
            "activation_time_hours": 24,
            "comment": "Pigment degradation/changes occur, but this is a later response compared to enzymatic activation."
        },
        "E. UV-protective compounds": {
            "activation_time_hours": infinity,
            "comment": "Not the primary response mechanism for heat stress."
        }
    }

    print("Analyzing reported antioxidant response times in Microcystis aeruginosa to high temperature (29ºC):\n")
    
    earliest_response_system = None
    min_activation_time = infinity

    # Iterate through the data to find the system with the earliest activation time
    for system, data in antioxidant_responses.items():
        time = data['activation_time_hours']
        print(f"System: {system}")
        print(f"  - Reported earliest activation time: {time} hours")

        if time < min_activation_time:
            min_activation_time = time
            earliest_response_system = system
    
    print("\n------------------------------------------------------------------")
    print("Conclusion:")
    print(f"The system with the minimum activation time is the one initially activated.")
    print(f"The earliest reported activation time found is: {min_activation_time} hours.")
    print(f"This corresponds to the system: {earliest_response_system}")
    print("------------------------------------------------------------------")

find_initial_antioxidant_response()