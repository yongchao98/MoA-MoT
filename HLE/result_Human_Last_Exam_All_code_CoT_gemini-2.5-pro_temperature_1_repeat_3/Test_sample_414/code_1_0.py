import numpy as np

# Step 1: Define the half-lives in days
T_p = 12.75  # Half-life of the parent, Ba-140
T_d = 1.678  # Half-life of the daughter, La-140

# Step 2: Calculate the decay constants (lambda) in day^-1
lambda_p = np.log(2) / T_p
lambda_d = np.log(2) / T_d

# Step 3: Define the activity ratio and time interval
activity_ratio = 2.1 / 1.4
time_interval = 14.0 # days

# Step 4: Solve the equation for T
# The equation derived from the ratio of activities is:
# activity_ratio = (exp(-lambda_p*(T+14)) - exp(-lambda_d*(T+14))) / (exp(-lambda_p*T) - exp(-lambda_d*T))
# This can be rearranged to:
# exp((lambda_d - lambda_p)*T) = (activity_ratio - exp(-lambda_d*14)) / (activity_ratio - exp(-lambda_p*14))

exp_term_p = np.exp(-lambda_p * time_interval)
exp_term_d = np.exp(-lambda_d * time_interval)

numerator = activity_ratio - exp_term_d
denominator = activity_ratio - exp_term_p

rhs = numerator / denominator
T = np.log(rhs) / (lambda_d - lambda_p)

# Step 5: Print the results and the equation steps
print("This problem tracks the grow-in of a daughter nuclide (La-140) from a parent (Ba-140).")
print("We assume the measured activity is dominated by the daughter, La-140.")
print("\nLet T be the time between separation and the first measurement.")
print("The ratio of the daughter's activities at T+14 days and T is equal to the ratio of measured activities.")
print(f"\nActivity Ratio = {2.1} / {1.4} = {activity_ratio}")
print("\nThe governing equation is:")
print("A(T+14)/A(T) = [exp(-λ_p*(T+14)) - exp(-λ_d*(T+14))] / [exp(-λ_p*T) - exp(-λ_d*T)] = 1.5")
print("\nSolving for T with the given half-lives:")
print(f"Half-life of Ba-140 (T_p): {T_p} days -> λ_p = {lambda_p:.4f} day^-1")
print(f"Half-life of La-140 (T_d): {T_d} days -> λ_d = {lambda_d:.4f} day^-1")
print("\nThe rearranged equation solved for T is:")
print("T = ln[(R - exp(-λ_d*14)) / (R - exp(-λ_p*14))] / (λ_d - λ_p)")
print(f"T = ln[({activity_ratio:.1f} - {exp_term_d:.4f}) / ({activity_ratio:.1f} - {exp_term_p:.4f})] / ({lambda_d:.4f} - {lambda_p:.4f})")
print(f"T = ln[{numerator:.4f} / {denominator:.4f}] / {lambda_d - lambda_p:.4f}")
print(f"T = ln[{rhs:.4f}] / {lambda_d - lambda_p:.4f}")
print(f"T = {np.log(rhs):.4f} / {lambda_d - lambda_p:.4f}")
print(f"\nThe calculated time T is: {T:.2f} days.")

print("\nTherefore, the approximate time between the sample irradiation (and subsequent separation) and the first analysis is about 1 day.")

# Final answer in the required format
# Round the final answer to one decimal place as it is an approximation
final_answer = round(T, 1)
# print(f"\n<<<{final_answer}>>>")