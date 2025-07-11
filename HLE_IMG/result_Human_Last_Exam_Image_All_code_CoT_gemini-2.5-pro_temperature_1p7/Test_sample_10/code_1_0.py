import math

# Step 1: Define the known data points.
# n is the number of 'acenaphtho' units fused to the core.
# B is the inversion barrier in kcal/mol.
n0, B0 = 0, 10.0
n2, B2 = 2, 49.0
n3 = 3

# Step 2: Formulate the model B(n) = a*n^2 + c and solve for the parameters.
# From the data point (n0, B0), we can find c.
# B(0) = a * 0^2 + c = B0 => c = B0
c = B0

# From the data point (n2, B2), we can find a.
# B(2) = a * n2^2 + c = B2 => a = (B2 - c) / n2^2
a = (B2 - c) / (n2**2)

# Step 3: Use the model to predict the barrier for n3.
# B(n3) = a * n3^2 + c
B3 = a * (n3**2) + c

# Step 4: Print the results and the final answer.
print("Predicting the inversion barrier based on the provided data.")
print("-" * 50)
print("Let 'n' be the number of fused acenaphtho units.")
print(f"Data point 1: For n = {n0}, the barrier is {B0} kcal/mol.")
print(f"Data point 2: For n = {n2}, the barrier is {B2} kcal/mol.")
print("-" * 50)
print("Assuming a quadratic model for the inversion barrier B(n) of the form: B(n) = a * n^2 + c")
print(f"Using n=0, B(0)=10: We find c = {c}")
print(f"Using n=2, B(2)=49: We find a = ({B2} - {c}) / {n2}^2 = {a}")
print(f"\nThe resulting model is: B(n) = {a} * n^2 + {c}")
print("-" * 50)
print(f"For triacenaphtho[...]triphenylene, n = {n3}.")
print("The predicted inversion barrier is:")
print(f"B({n3}) = {a} * {n3}^2 + {c} = {B3} kcal/mol")
print("-" * 50)
print("As requested, the value with an accuracy of integers is:")
print(f"Predicted inversion barrier = {math.ceil(B3)} kcal/mol")

# Outputting the final answer in the required format.
# Note: The calculation yields 97.75, which rounds to 98. However, since the question asks for accuracy of integers and chemists often report barriers as integers or with one decimal place, rounding up (ceil) is a reasonable choice given the increasing trend. Rounding to nearest would also give 98.
final_answer = math.ceil(B3)
# final_answer_str = f"<<<{final_answer}>>>"
# This final print statement is for the platform, the value is calculated above.