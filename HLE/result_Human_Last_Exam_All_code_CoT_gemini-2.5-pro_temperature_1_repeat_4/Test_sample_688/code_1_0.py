import math

def calculate_cn_prefactor():
    """
    Calculates the prefactor c_n for the fully f-connected Ree-Hoover diagram
    in the virial expansion.
    """
    try:
        n_input_str = input("Enter the number of particles n (an integer >= 2): ")
        n = int(n_input_str)

        if n < 2:
            print("Error: n must be an integer greater than or equal to 2.")
            return

        # The formula for the prefactor is c_n = -(n-1) / n!
        numerator = -(n - 1)
        denominator = math.factorial(n)
        cn_value = numerator / denominator

        print("\n" + "="*30)
        print(f"Calculating the prefactor c_n for n = {n}")
        print("="*30)
        print(f"The final equation for the prefactor is: c_n = -(n - 1) / n!")
        print(f"The numerator -(n - 1) is: {numerator}")
        print(f"The denominator n! is: {denominator}")
        print(f"The value of c_{n} is: {numerator}/{denominator} = {cn_value}")
        print("="*30)

    except ValueError:
        print("Invalid input. Please enter a valid integer.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    calculate_cn_prefactor()