def find_initial_antioxidant_response():
    """
    This function models biological knowledge to determine the initial
    antioxidant system activated in response to acute stress like high temperature.
    """

    # A simplified knowledge base describing antioxidant systems in cyanobacteria.
    # The 'response_speed' key is crucial for identifying the 'initial' response.
    antioxidant_systems = {
        "A": {
            "name": "Liposoluble antioxidants",
            "description": "Includes compounds like α-tocopherol (Vitamin E) and carotenoids. They are located in membranes and lipid bodies.",
            "role": "Scavenge radicals within membranes, protect lipids from peroxidation.",
            "response_speed": "constitutive/secondary"
        },
        "B": {
            "name": "Hydrosoluble antioxidants",
            "description": "Includes compounds like ascorbic acid (Vitamin C) and glutathione.",
            "role": "Scavenge radicals in the cytoplasm and regenerate other antioxidants.",
            "response_speed": "secondary/regenerative"
        },
        "C": {
            "name": "Enzymatic antioxidants",
            "description": "Includes enzymes like Superoxide Dismutase (SOD), Catalase (CAT), and Ascorbate Peroxidase (APX).",
            "role": "Primary enzymatic detoxification of Reactive Oxygen Species (ROS). SOD converts superoxide to H2O2, and CAT/APX convert H2O2 to water.",
            "response_speed": "primary/rapid"
        },
        "D": {
            "name": "Photosynthetic pigments",
            "description": "Includes chlorophylls and phycobiliproteins.",
            "role": "While they can quench some ROS, they are also a primary source of ROS under stress and are often damaged (photo-oxidation). Their primary role is light harvesting.",
            "response_speed": "passive/secondary"
        },
        "E": {
            "name": "UV-protective compounds",
            "description": "Includes mycosporine-like amino acids (MAAs) and scytonemin.",
            "role": "Primarily screen harmful UV radiation, not a direct response to thermal stress.",
            "response_speed": "slow/adaptive"
        }
    }

    # The stress condition is high temperature, which causes acute oxidative stress.
    # We are looking for the system with the most rapid, primary response.
    initial_responder = None
    for key, system in antioxidant_systems.items():
        if "primary/rapid" in system["response_speed"]:
            initial_responder = system
            break

    if initial_responder:
        print("Analysis of Antioxidant Systems:")
        print("-" * 30)
        for key, system in antioxidant_systems.items():
            print(f"System ({key}): {system['name']}")
            print(f"  Role: {system['role']}")
            print(f"  Response Type: {system['response_speed']}")
            print("-" * 30)

        print("\nConclusion:")
        print(f"The initial response to acute oxidative stress from high temperature is managed by the system designed for rapid, primary detoxification of reactive oxygen species.")
        print(f"This corresponds to the '{initial_responder['name']}'.")
        print("\nStudies on Microcystis aeruginosa confirm that the activity of enzymes like Superoxide Dismutase (SOD) and Catalase (CAT) increases rapidly upon high temperature exposure to counteract the initial burst of oxidative stress.")

    else:
        print("Could not determine the initial responder based on the defined knowledge base.")

find_initial_antioxidant_response()