# The established order of enclitics in Old Russian follows a specific hierarchy.
# The given enclitics are: бы, же, еси, бо, мя.

# Let's classify them and arrange them according to the linguistic rules:
# 1. 'бо' (conjunction, "for, because")
# 2. 'мя' (accusative personal pronoun, "me")
# 3. 'еси' (2nd person singular of the verb 'to be')
# 4. 'бы' (conditional particle)
# 5. 'же' (emphatic particle, always last in the chain)

correct_order = ["бо", "мя", "еси", "бы", "же"]

# Let's construct a hypothetical sentence fragment to illustrate:
# "Рече бо мя еси бы же..." meaning something like "For you would have told me..."
# The script will now print the final, correct sequence.

print("The correct order of the enclitics is:")
# We use ' -> ' to represent the sequence they would appear in.
final_sequence = " -> ".join(correct_order)
print(final_sequence)