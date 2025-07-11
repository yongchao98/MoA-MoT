import textwrap

def find_primary_condition():
    """
    Analyzes the link between Doppelgänger syndrome and neurological conditions
    to identify the primary disease identified by modern neuroscience.
    """
    question = "With neuroscientific development and improvements in imaging diagnosis, a specific condition was identified as the primary disease underlying the doppelgänger syndrome. What is it?"

    # A knowledge base representing conditions that can cause autoscopic phenomena.
    # The 'primary' flag is based on neuroscientific findings linking specific,
    # observable brain events (like seizures in the temporoparietal junction) to the syndrome.
    knowledge_base = {
        "Schizophrenia": {"primary": False, "notes": "Can cause a wide range of hallucinations, but is a psychiatric rather than a focal neurological cause."},
        "Migraine with aura": {"primary": False, "notes": "Can cause complex visual disturbances, but is less commonly the direct cause of a full autoscopic experience."},
        "Stroke or Tumor": {"primary": False, "notes": "Lesions in specific brain areas (like the temporoparietal junction) can cause the syndrome, but they are the structural cause, not the disease process itself in all cases."},
        "Epilepsy": {"primary": True, "notes": "Focal seizures, particularly in the temporoparietal or occipital lobes, are a classic and well-documented cause. Modern imaging (EEG, fMRI) can pinpoint this epileptic activity as the direct trigger for the experience."}
    }

    print(f"Query: {question}\n")
    print("Searching neuroscientific knowledge base...")

    primary_condition = None
    for condition, data in knowledge_base.items():
        if data["primary"]:
            primary_condition = condition
            notes = data["notes"]
            break

    if primary_condition:
        print(f"\nPrimary condition identified: {primary_condition}")
        print("\nExplanation:")
        wrapped_notes = textwrap.fill(notes, width=80)
        print(wrapped_notes)
    else:
        print("Could not identify a primary condition based on the available data.")

find_primary_condition()
