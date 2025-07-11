def find_doppelganger_condition():
    """
    Simulates searching a knowledge base to find the neurological syndrome
    associated with the doppelgänger phenomenon.
    """
    # A simple knowledge base mapping conditions to descriptive keywords.
    knowledge_base = {
        "Heautoscopy": ["seeing own body", "reflection", "symptom"],
        "Fregoli Delusion": ["disguise", "single person", "syndrome"],
        "Temporal Lobe Epilepsy": ["seizures", "hallucinations", "cause"],
        "Capgras Syndrome": ["impostor", "double", "replaced", "syndrome"]
    }

    # Keywords derived from the query for a "doppelgänger" (double) that is a "primary disease" (syndrome).
    query_keywords = ["double", "syndrome"]

    found_condition = "Condition not found in knowledge base."

    # Iterate through the knowledge base to find the best match.
    for condition, keywords in knowledge_base.items():
        # Check if the condition's keywords contain all the query keywords.
        if all(query_keyword in keywords for query_keyword in query_keywords):
            found_condition = condition
            break  # Stop after finding the first match

    print("Searching for a primary condition (syndrome) related to seeing a 'double'...")
    print(f"The identified condition is: {found_condition}")

find_doppelganger_condition()