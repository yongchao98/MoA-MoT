def sort_old_russian_enclitics():
    """
    Sorts a given list of Old Russian enclitics based on their established grammatical order.
    """
    # The established order of enclitics (or a representative sample).
    # 1. Conditional particle (бы)
    # 2. Pronominal accusative/dative (мя)
    # 3. Verb "to be" (еси)
    # 4. Other particles/conjunctions (же, бо). Within this group, 'же' precedes 'бо'.
    correct_order = ["бы", "мя", "еси", "же", "бо"]

    # The list of enclitics provided by the user.
    enclitics_to_sort = ["бы", "же", "еси", "бо", "мя"]

    # Sort the user's list based on the index in the correct_order list.
    sorted_enclitics = sorted(enclitics_to_sort, key=lambda x: correct_order.index(x))

    print("The correct order for the given enclitics is:")
    # The final print statement will output each item in the ordered sequence.
    print(" -> ".join(sorted_enclitics))

sort_old_russian_enclitics()

# The final answer in the requested format.
# Each item from the final ordered list is presented.
print("\n<<<бы, мя, еси, же, бо>>>")