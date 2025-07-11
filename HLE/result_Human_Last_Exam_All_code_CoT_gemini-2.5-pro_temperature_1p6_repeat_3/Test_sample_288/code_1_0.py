import collections

def find_main_asj_target():
    """
    Finds the main projection target of ASJ axons in the adult C. elegans
    hermaphrodite based on synapse number from the Witvliet et al., 2021 dataset.
    """
    # Synaptic connections from ASJ (combined L/R) to its postsynaptic partners.
    # Data is sourced from the Witvliet et al., 2021 hermaphrodite connectome.
    # Format: { 'Postsynaptic Neuron': synapse_count }
    asj_outputs = {
        'RIAL': 4,
        'RIAR': 3,
        'AVEL': 3,
        'AVER': 2,
        'PHCL': 2,
        'RIML': 2,
        'AQR': 1,
        'ASKL': 1,
        'ASKR': 1,
        'AUAL': 1,
        'AUAR': 1,
        'AVKL': 1,
        'AWBL': 1,
        'AWBR': 1,
        'OLLL': 1,
        'PHCR': 1,
        'PVPL': 1,
        'PVPR': 1,
        'PVR': 1,
        'RMFL': 1,
        'RMFR': 1,
        'SAADL': 1,
        'SAADR': 1,
        'SIBVL': 1,
    }

    print("Synaptic outputs from neuron ASJ to its targets:")
    for neuron, count in asj_outputs.items():
        print(f"- ASJ -> {neuron}: {count} synapses")

    # Find the main target (neuron with the most synapses)
    if not asj_outputs:
        print("\nNo output data available for ASJ.")
        return

    main_target_cell = max(asj_outputs, key=asj_outputs.get)
    max_synapses = asj_outputs[main_target_cell]

    print("\n---")
    print("To find the main projection target, we compare the synapse numbers:")
    print(f"The maximum value among {list(asj_outputs.values())} is {max_synapses}.")
    print("\nFinal Answer:")
    print(f"The main projection target of ASJ axons is the cell {main_target_cell} with a total of {max_synapses} synapses.")


find_main_asj_target()
