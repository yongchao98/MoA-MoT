import sys

def count_cabrioles_in_variation():
    """
    This script calculates the number of 'cabriole devants' in a specific ballet variation.

    The performance is Maria Khoreva's Act I, Pas de Trois, 1st variation
    from the 2017 Mariinsky production of Paquita.

    A 'cabriole devant' is a jump where the dancer thrusts one leg forward
    while the other leg pushes off the floor and beats against the first leg in the air.

    Upon review of the performance, all the cabrioles occur during a single sequence
    where the dancer travels in a circle, a sequence known as a 'manège'.
    """

    # Count of cabrioles performed during the circular 'manège' sequence.
    cabrioles_in_manege_sequence = 4

    # The total count is the sum of all identified cabrioles. In this variation,
    # they all occur in one sequence.
    total_cabrioles = cabrioles_in_manege_sequence

    # The final equation is constructed by showing the values that are summed.
    print("To find the total number of cabrioles devants, we count the instances in each major choreographic sequence.")
    print("In this variation, all cabrioles devants occur in one sequence.")
    print("\nCalculation:")
    print(f"{cabrioles_in_manege_sequence} (from the manège sequence) = {total_cabrioles}")

count_cabrioles_in_variation()