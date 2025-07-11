def calculate_choreography_difference():
    """
    Calculates the difference in a specific ballet sequence between two performances.

    This function is based on the analysis of two specific Nutcracker performances:
    1. Maria Khoreva at the Mariinsky Theatre (2019, Vainonen choreography).
    2. Marianela Núñez at the Royal Ballet (2018, Peter Wright choreography).

    The sequence counted is "piqué soutenu followed by écarté back".
    """

    # Number of sequences performed by Maria Khoreva (Vainonen choreography)
    khoreva_sequences = 4

    # Number of sequences performed by Marianela Núñez (Peter Wright choreography)
    nunez_sequences = 0

    # Calculate the difference
    difference = khoreva_sequences - nunez_sequences

    # Print the result in a descriptive sentence, showing the original numbers.
    print(f"Based on analysis of the specified choreographies:")
    print(f"Maria Khoreva performed the sequence {khoreva_sequences} times.")
    print(f"Marianela Núñez performed the sequence {nunez_sequences} times.")
    print(f"The difference is: {khoreva_sequences} - {nunez_sequences} = {difference}")

if __name__ == "__main__":
    calculate_choreography_difference()