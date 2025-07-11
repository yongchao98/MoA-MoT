import math

def calculate_mistake_bound():
    """
    Calculates and explains the upper bound on mistakes for a variant of the experts problem.
    """
    try:
        n_str = input("Enter the total number of experts (n): ")
        c_str = input("Enter the mistake threshold for removing an expert (c): ")

        n = int(n_str)
        c = int(c_str)

        if n <= 1:
            print("Error: Number of experts (n) must be greater than 1.")
            return
        if c <= 1:
            print("Error: Mistake threshold (c) must be greater than 1.")
            return

        # Bound for M1: Mistakes where the true expert is also wrong.
        # The true expert makes at most c-1 mistakes.
        m1_bound = c - 1

        # Bound for M2: Mistakes where the algorithm is wrong but the true expert is right.
        # Total mistakes available to n-1 false experts is (n-1)*c.
        # Each M2 mistake requires at least 2 false expert mistakes.
        total_false_mistakes = (n - 1) * c
        m2_bound = total_false_mistakes / 2

        # The total bound is the sum of the two parts.
        total_bound = m1_bound + m2_bound

        print("\n--- Upper Bound Calculation ---")
        print(f"The final formula for the upper bound on mistakes M is: M <= (c - 1) + ((n - 1) * c) / 2")

        print("\nStep 1: Calculate the bound for mistakes where the true expert is also wrong (M1).")
        print(f"M1 <= c - 1")
        print(f"M1 <= {c} - 1 = {m1_bound}")

        print("\nStep 2: Calculate the bound for mistakes where the true expert is right (M2).")
        print(f"M2 <= ((n - 1) * c) / 2")
        print(f"M2 <= (({n} - 1) * {c}) / 2 = ({n-1} * {c}) / 2 = {total_false_mistakes} / 2 = {m2_bound}")

        print("\nStep 3: Calculate the total upper bound by summing the bounds for M1 and M2.")
        print(f"M <= {m1_bound} + {m2_bound}")
        
        # Note: The number of mistakes must be an integer. The bound can be a non-integer.
        # However, often this is expressed as floor(bound). Let's present the precise value.
        print(f"\nFinal Upper Bound: {total_bound}")


    except ValueError:
        print("Invalid input. Please enter integer values for n and c.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    calculate_mistake_bound()