import collections

# Step 1: Identify the main choreographic phrase containing the cabrioles.
# In this variation, Maria Khoreva performs a sequence of four gargouillades, 
# each followed by a cabriole devant, as she travels down a diagonal.

# Step 2: Count the number of cabrioles in each repetition of the phrase.
cabriole_counts = [
    1,  # First phrase: gargouillade, cabriole
    1,  # Second phrase: gargouillade, cabriole
    1,  # Third phrase: gargouillade, cabriole
    1   # Fourth phrase: gargouillade, cabriole
]

# Step 3: Calculate the total sum.
total_cabrioles = sum(cabriole_counts)

# Step 4: Construct and print the equation and the result.
# The following code will display each number in the final equation.
equation_parts = [str(count) for count in cabriole_counts]
equation_string = " + ".join(equation_parts)
print(f"The total number of cabrioles devants is calculated as follows:")
print(f"{equation_string} = {total_cabrioles}")
