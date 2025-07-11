from fractions import Fraction

def predict_inversion_barrier():
    """
    Calculates the inversion barrier for a PAH based on a linear model
    derived from two other molecules in the same series.
    """
    # Data for the first molecule (corannulene)
    n1 = 20  # Number of carbon atoms
    b1 = 10  # Inversion barrier in kcal/mol

    # Data for the second molecule (diacenaphthochrysene derivative)
    n2 = 38  # Number of carbon atoms
    b2 = 49  # Inversion barrier in kcal/mol

    # Data for the third molecule (triacenaphthotriphenylene derivative)
    n3 = 42  # Number of carbon atoms

    print("Step 1: Establishing a linear model B = a*N + b for the inversion barrier (B) based on the number of carbon atoms (N).")
    print(f"Molecule 1: N1 = {n1}, B1 = {b1} kcal/mol")
    print(f"Molecule 2: N2 = {n2}, B2 = {b2} kcal/mol")
    print("-" * 30)

    # Calculate the slope 'a' of the linear model
    # To maintain precision, we use the Fraction class
    a = Fraction(b2 - b1, n2 - n1)

    # Calculate the intercept 'b' of the linear model
    b = Fraction(b1) - a * Fraction(n1)

    print("Step 2: Calculating the parameters of the linear model.")
    print(f"Slope a = (B2 - B1) / (N2 - N1) = ({b2} - {b1}) / ({n2} - {n1}) = {b2 - b1}/{n2 - n1} = {a}")
    print(f"Intercept b = B1 - a * N1 = {b1} - ({a}) * {n1} = {b}")
    print(f"The resulting linear equation is: B = ({a})*N + ({b})")
    print("-" * 30)

    # Predict the barrier for the third molecule
    b3 = a * Fraction(n3) + b
    b3_rounded = round(float(b3))

    print(f"Step 3: Predicting the inversion barrier for Molecule 3 with N3 = {n3} carbon atoms.")
    print(f"Calculation: B3 = ({a}) * {n3} + ({b})")
    print(f"B3 = {float(a * Fraction(n3))} + ({float(b)})")
    print(f"B3 = {float(b3):.2f} kcal/mol")
    print("-" * 30)

    print("Final Answer:")
    print(f"The predicted inversion barrier for triacenaphtho[3,2,1,8-cdefg:3',2',1',8'-ijklm:3'',2'',1'',8''-opqra]triphenylene is {b3_rounded} kcal/mol (rounded to the nearest integer).")

if __name__ == '__main__':
    predict_inversion_barrier()
<<<58>>>