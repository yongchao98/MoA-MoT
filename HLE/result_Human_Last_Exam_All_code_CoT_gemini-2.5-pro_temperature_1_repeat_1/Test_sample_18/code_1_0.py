import numpy as np

# The problem states that the K-matrix for the bosonic integer quantum Hall
# state at nu=2 is the Pauli matrix sigma_x.
K_boson = np.array([[0, 1], 
                    [1, 0]])

# The problem also states that these bosons are formed from fermions by attaching
# two flux quanta to each fermion. For a two-component system, this corresponds
# to a flux attachment matrix M.
M = np.array([[2, 0], 
              [0, 2]])

# The K-matrix of the original fermionic state, K_fermion, is related to the
# K-matrix of the composite bosons, K_boson, by the formula:
# K_boson = K_fermion + M
# Therefore, to find K_fermion, we perform a matrix subtraction:
# K_fermion = K_boson - M
K_fermion = K_boson - M

# Print the final equation, showing each number in the calculation as requested.
print("The K-matrix of the resulting fractional state is K_fermion.")
print("The calculation is: K_fermion = K_boson - M")
print("")
print("K_fermion =")
# Print the matrix equation element by element for clarity
# The f-string formatting helps align the numbers.
print(f"| {K_boson[0,0]:^3} {K_boson[0,1]:^3} | - | {M[0,0]:^3} {M[0,1]:^3} | = | {K_fermion[0,0]:^3} {K_fermion[0,1]:^3} |")
print(f"| {K_boson[1,0]:^3} {K_boson[1,1]:^3} |   | {M[1,0]:^3} {M[1,1]:^3} |   | {K_fermion[1,0]:^3} {K_fermion[1,1]:^3} |")