import sys
import io

# Helper function to redirect stdout to capture print statements
def capture_output(func, *args, **kwargs):
    old_stdout = sys.stdout
    sys.stdout = new_stdout = io.StringIO()
    try:
        func(*args, **kwargs)
        return new_stdout.getvalue()
    finally:
        sys.stdout = old_stdout

def calculate_selection_coefficient():
    """
    This function demonstrates how to calculate the cost of gene flow (selection coefficient) in yeast.
    It uses placeholder values for fitness measurements.
    """

    # Step 1: Define fitness values from a hypothetical experiment.
    # Fitness can be measured by growth rate, biomass, etc.
    # We assume the parental line is the fitter reference genotype.
    fitness_parent = 1.0  # Fitness of the non-hybridized parental line (control).
    fitness_hybrid = 0.9  # Measured fitness of the hybrid line (a 10% reduction).

    # Step 2: Calculate the relative fitness (W) of the hybrid.
    # This compares the hybrid's fitness to the parent's.
    relative_fitness_W = fitness_hybrid / fitness_parent

    # Step 3: Calculate the selection coefficient (s), which represents the "cost" of gene flow.
    selection_coefficient_s = 1 - relative_fitness_W

    # Step 4: Print the entire process, including the final equation with numbers.
    print("The cost of gene flow is measured by the selection coefficient (s).")
    print("The formula is s = 1 - W, where W is the relative fitness of the hybrid.\n")

    print("--- Calculation Steps ---")
    print(f"1. Measure fitness of parental line (no gene flow): {fitness_parent}")
    print(f"2. Measure fitness of hybrid line (with gene flow): {fitness_hybrid}\n")

    print(f"3. Calculate relative fitness (W) of the hybrid:")
    print(f"   W = fitness_hybrid / fitness_parent = {fitness_hybrid} / {fitness_parent} = {relative_fitness_W}\n")
    
    print("4. Calculate the selection coefficient (s) using the final equation:")
    # This prints the final equation with each number explicitly shown.
    print(f"   s = 1 - W")
    print(f"   s = 1 - {relative_fitness_W}")
    print(f"   s = {selection_coefficient_s}\n")
    
    print(f"The resulting cost due to gene flow (selection coefficient) is {selection_coefficient_s:.2f}.")

# Execute the function and print its output to the console.
output = capture_output(calculate_selection_coefficient)
print(output)