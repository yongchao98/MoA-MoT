# Step 1: Evaluate the multiplication inside the first parentheses
step1 = -483 * -873

# Step 2: Evaluate the division inside the second parentheses
step2 = -802 / 2

# Step 3: Evaluate the addition inside the second parentheses
step3 = -253 + step2

# Step 4: Evaluate the negation of the results from steps 1 and 3
step4a = -step1
step4b = -step3

# Step 5: Evaluate the multiplication of the results from step 4
step5 = step4a * step4b

# Step 6: Add 383 to the result from step 5
result = step5 + 383

print(result)