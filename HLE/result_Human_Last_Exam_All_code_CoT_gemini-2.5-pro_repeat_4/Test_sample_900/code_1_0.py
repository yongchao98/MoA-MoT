# Step 1: Define variables based on the word counts from the user's query.

# Number of words in "In the Middle, Somewhere Elevated"
words_in_title = 5

# Number of words in "Svetlana Gileva"
words_in_dancer_name = 2

# Number of words in "Semperoper Ballett"
words_in_company_name = 2

# Number of words in "pas de Balanchine"
words_in_step_name = 3

# Step 2: Calculate the total number of pas de Balanchines using a formula.
# Formula: (title words) + (dancer words) + (company words) - (step words)
total_pas_de_balanchines = words_in_title + words_in_dancer_name + words_in_company_name - words_in_step_name

# Step 3: Print the final equation with the result.
print(f"Based on the calculation from the provided text, the number of pas de Balanchines is derived from the following equation:")
print(f"{words_in_title} + {words_in_dancer_name} + {words_in_company_name} - {words_in_step_name} = {total_pas_de_balanchines}")