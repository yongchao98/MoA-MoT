def analyze_sorting_efficiency():
    """
    Analyzes the efficiency of a cell sorting experiment and determines the most likely cause for low yield.
    """
    total_sorted_wells = 95
    successful_wells = 4

    # Calculate the efficiency
    efficiency = (successful_wells / total_sorted_wells) * 100

    # Print the step-by-step calculation
    print("Experimental Efficiency Calculation:")
    print(f"Total sorted wells (potential doublets): {total_sorted_wells}")
    print(f"Wells with paired heavy and light chains: {successful_wells}")
    print(f"Efficiency = ({successful_wells} / {total_sorted_wells}) * 100 = {efficiency:.2f}%")
    print("\n")

    # Explain the reasoning for the low efficiency
    print("Analysis of Low Efficiency:")
    print("The calculated efficiency of ~4.21% is extremely low. Let's evaluate the most likely cause:")
    print(
        "The problem states that double-positive (RFP+ FITC+) events were sorted. However, a major challenge in flow cytometry is distinguishing between true doublets (one tumor cell and one B cell physically bound) and coincidence events (one tumor cell and one B cell passing through the laser beam in very close succession but not physically attached)."
    )
    print(
        "A sorter cannot easily distinguish between these two scenarios and will sort both as 'double-positive' events. If true, stable interactions are rare, the vast majority of events gated as 'doublets' will actually be these coincidence events."
    )
    print(
        "When such a 'false doublet' is sorted, the well may receive only a single tumor cell, a single B cell, or one of each that weren't truly interacting. If a well receives only a tumor cell, the subsequent PCR for B cell heavy and light chains will fail, yielding no band. This scenario perfectly explains why 91 out of 95 wells produced nothing."
    )
    print(
        "While other factors like RNA competition (E) could reduce efficiency, they are unlikely to cause a >95% failure rate, especially with a sensitive nested PCR protocol. The most direct and common explanation for this specific outcome is a failure at the cell isolation step."
    )

analyze_sorting_efficiency()
print("<<<D>>>")