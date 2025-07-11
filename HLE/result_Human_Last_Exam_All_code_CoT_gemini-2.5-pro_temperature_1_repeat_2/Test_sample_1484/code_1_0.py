def calculate_ballet_sequence_difference():
    """
    Calculates the difference in a specific ballet sequence between two performances.
    """

    # Count for Marianela Núñez at the Royal Ballet (2018, Peter Wright choreography)
    # In this version, the sequence is performed four times in a row.
    nunez_count = 4

    # Count for Maria Khoreva at the Mariinsky Theatre (2019, Vasily Vainonen choreography)
    # This specific sequence is not present in the Vainonen choreography.
    khoreva_count = 0

    # Calculate the difference
    difference = nunez_count - khoreva_count

    # Print the findings and the final calculation
    print(f"In the 2018 Royal Ballet performance, Marianela Núñez performed the sequence {nunez_count} times.")
    print(f"In the 2019 Mariinsky Theatre performance, Maria Khoreva performed the sequence {khoreva_count} times.")
    print(f"The difference is {nunez_count} - {khoreva_count} = {difference}.")

calculate_ballet_sequence_difference()
<<<4>>>