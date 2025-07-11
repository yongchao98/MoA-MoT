def analyze_experiment_control():
    """
    Analyzes the experimental design and determines the role and timing of the anti-FLAG antibody.
    """
    
    # Experimental parameters
    inhibitor = "GalNAc"
    inhibitor_concentration_mM = 500
    protein_tag = "FLAG tag"
    primary_antibody = "Anti-MUC1/TN-antigen antibody"
    control_antibody = "Anti-FLAG antibody"

    print("### Experimental Analysis ###")
    print(f"The goal is to show that an antibody specifically binds a glycan on MUC1.")
    print(f"A competitive inhibition experiment is used with {inhibitor_concentration_mM} mM of free {inhibitor}.")
    
    print("\n### The Problem with the Inhibitor ###")
    print(f"A high concentration like {inhibitor_concentration_mM} mM could have off-target effects.")
    print("A major concern is that it could cause the cell to reduce the amount of MUC1 protein on its surface.")
    print("This would lead to a false positive result for competitive inhibition.")

    print(f"\n### The Role of the {control_antibody} ###")
    print(f"The MUC1 protein has a {protein_tag}. The {control_antibody} binds to this tag, not the sugar.")
    print("Therefore, its binding is independent of the inhibitor.")
    print(f"By measuring the {control_antibody} signal, we can confirm that the total amount of MUC1 protein on the cell surface is not changed by the {inhibitor_concentration_mM} mM {inhibitor} treatment.")

    print("\n### Timing of Antibody Addition ###")
    print(f"The {control_antibody} binds directly to the target protein on the cell surface.")
    print("In immunolabeling protocols, antibodies that bind directly to the antigen are called 'primary antibodies'.")
    print("Therefore, it must be added during the primary antibody incubation step.")

    print("\n### Conclusion ###")
    print("The anti-flag antibody is added with the primary antibodies to verify that the high concentration of GalNAc has not altered the surface expression of the MUC1 protein itself.")
    print("This directly corresponds to Answer Choice C.")

# Run the analysis
analyze_experiment_control()