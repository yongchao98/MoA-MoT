# Define the truth values for each person
fidel = False  # Fidel lies
fletcher = not fidel  # Fletcher says Fidel tells the truth
yoland = not fletcher  # Yoland says Fletcher lies
raymond = not yoland  # Raymond says Yoland lies
leda = not raymond  # Leda says Raymond lies

# Print whether Leda tells the truth
print(leda)