def interpret_histogram(antibody_name, has_overlap):
    """
    This function interprets a simplified model of a flow cytometry histogram.

    Args:
        antibody_name (str): The name of the antibody ('A', 'B', or 'C').
        has_overlap (bool): True if the stained population overlaps with the unstained control.

    Returns:
        bool: True if all cells were stained, False otherwise.
    """
    print(f"--- Analyzing Antibody {antibody_name} ---")
    if has_overlap:
        print("Result: The stained (red) and unstained (black) populations overlap.")
        print(f"This means not all cells were stained by antibody {antibody_name}.")
        return False
    else:
        print("Result: The stained (red) population is completely separate from the unstained (black) population.")
        print(f"This means all cells were stained by antibody {antibody_name}.")
        return True

def main():
    """
    Based on visual inspection, we determine if there is an overlap for each experiment.
    """
    print("Analyzing which antibody stained all cells in the sample...\n")

    # From visual analysis of Histogram A, there is an overlap.
    a_stained_all = interpret_histogram('A', has_overlap=True)
    print("")

    # From visual analysis of Histogram B, there is no overlap.
    b_stained_all = interpret_histogram('B', has_overlap=False)
    print("")

    # From visual analysis of Histogram C, there is an overlap.
    c_stained_all = interpret_histogram('C', has_overlap=True)
    print("\n--- Summary ---")

    final_result = []
    if a_stained_all:
        final_result.append('A')
    if b_stained_all:
        final_result.append('B')
    if c_stained_all:
        final_result.append('C')
    
    if len(final_result) > 0:
        print(f"The antibody that stained all cells is: {', '.join(final_result)}")
    else:
        print("None of the antibodies stained all the cells.")


if __name__ == "__main__":
    main()