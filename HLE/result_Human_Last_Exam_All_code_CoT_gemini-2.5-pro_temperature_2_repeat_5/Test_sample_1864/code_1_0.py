def calculate_flow_controls():
    """
    Calculates the number of essential controls for a flow cytometry experiment.
    """
    # Number of channels/fluorophores in the experiment
    num_channels = 5

    # You always need one unstained control to measure autofluorescence.
    unstained_control = 1

    # You need one single-color compensation control for each channel.
    compensation_controls = num_channels

    # The total number of essential controls is the sum of the two types.
    total_controls = unstained_control + compensation_controls

    print("To determine the number of essential controls, we sum:")
    print("1. The unstained control (for baseline autofluorescence).")
    print("2. The single-color controls (one for each of the {} channels for compensation).".format(num_channels))
    print("\nCalculation:")
    print(f"{unstained_control} + {compensation_controls} = {total_controls}")
    print("\nTherefore, you should prepare a total of {} essential controls.".format(total_controls))

if __name__ == '__main__':
    calculate_flow_controls()