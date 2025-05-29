from sympy import symbols, Eq, solve

# Define the variable
QS = symbols('QS')

# Given values
RS = 9
cos_R = 4/9

# QR in terms of QS
QR = cos_R * QS

# Pythagorean theorem
equation = Eq(QS**2, QR**2 + RS**2)

# Solve for QS
QS_value = solve(equation, QS)
print(QS_value)